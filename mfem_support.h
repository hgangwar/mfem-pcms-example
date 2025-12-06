#ifndef mfem_support_h
#define mfem_support_h

#include "mfem.hpp"

namespace mfem_support
{
struct FEMSystem
{
  mfem::Mesh* mesh = nullptr;
  mfem::ParMesh* pmesh = nullptr;
  mfem::H1_FECollection* fec = nullptr;
  mfem::ParFiniteElementSpace* fes = nullptr;
  mfem::Array<int> ess_bdr;
  mfem::ParBilinearForm* a = nullptr;
  mfem::ParLinearForm* b = nullptr;
  mfem::ParGridFunction* x = nullptr;
  mfem::DomainLFIntegrator* rhs_int = nullptr;
};
// ------------------------------------------------------------
// Assign element attributes based on x-coordinate:
//   attribute = 1 if x ≤ 0.4
//   attribute = 2 if 0.4 < x <= 0.6
//   attribute = 3 otherwise
// ------------------------------------------------------------
void AssignAttributesByX(mfem::ParMesh& pmesh)
{
  for (int e = 0; e < pmesh.GetNE(); e++) {
    mfem::Vector center;
    pmesh.GetElementCenter(e, center);
    const double x = center[0];

    if (x <= 0.4)
      pmesh.SetAttribute(e, 1); // region 1 client_A only
    else if (x <= 0.6 && x > 0.4)
      pmesh.SetAttribute(e, 2); // region 2 overlap
    else
      pmesh.SetAttribute(e, 3); // region 3 client_B only
  }

  // Make sure ghost layer has the same attributes
  pmesh.ExchangeFaceNbrData();
}

// ------------------------------------------------------------
// Select DOFs that lie on the interior plane x = 0.5
// Produces an essential TRUE dof list: ess_tdofs
// ------------------------------------------------------------
void MarkInteriorPlaneDOFs(const mfem::ParMesh& pmesh,
                           const mfem::ParFiniteElementSpace& pfes,
                           mfem::Array<int>& ess_tdofs)
{
  const double xc = 0.5;
  const double tol = 1e-12;

  mfem::Array<int> vdof_marker(pfes.GetVSize());
  vdof_marker = 0;

  for (int v = 0; v < pmesh.GetNV(); v++) {
    const double* vx = pmesh.GetVertex(v);
    const double x = vx[0];

    if (std::fabs(x - xc) < tol) {
      mfem::Array<int> vdofs;
      pfes.GetVertexDofs(v, vdofs);

      for (int j = 0; j < vdofs.Size(); j++)
        vdof_marker[vdofs[j]] = 1;
    }
  }

  // convert local VDOF markers → TRUE DOFs in parallel
  pfes.GetEssentialTrueDofs(vdof_marker, ess_tdofs);
}

FEMSystem Init_FEMSystem(MPI_Comm comm, const std::string& mesh_file,
                         int order, const char client)
{
    FEMSystem sys;

    // Load and parallelize mesh
    sys.mesh  = new mfem::Mesh(mesh_file.c_str(), 1, 1, true);
    sys.pmesh = new mfem::ParMesh(comm, *sys.mesh);
    int dim   = sys.pmesh->Dimension();

    // Assign element attributes (1,2,3 based on x location)
    AssignAttributesByX(*sys.pmesh);

    // Attribute mask based on client
    const int max_attr = sys.pmesh->attributes.Max();
    mfem::Array<int> attr_mask(max_attr);
    attr_mask = 0;

    switch (client) {
        case 'A': attr_mask[0] = 1; attr_mask[1] = 1; break;  // 1 & 2
        case 'B': attr_mask[1] = 1; attr_mask[2] = 1; break;  // 2 & 3
        case 'M': attr_mask = 1; break;                       // 1,2,3
        default: throw std::invalid_argument("Unknown client type");
    }

    // -----------------------
    // Boundary conditions
    // -----------------------
    sys.ess_bdr.SetSize(sys.pmesh->bdr_attributes.Max());
    sys.ess_bdr = 0;

    int left_attr  = 1;  // must match Mesh bdr attributes
    int right_attr = 2;

    sys.ess_bdr[left_attr - 1]  = 1;
    sys.ess_bdr[right_attr - 1] = 1;

    // FE space
    sys.fec = new mfem::H1_FECollection(order, dim);
    sys.fes = new mfem::ParFiniteElementSpace(sys.pmesh, sys.fec);

    // ------------------------
    // Bilinear form (diffusion)
    // ------------------------
    sys.a = new mfem::ParBilinearForm(sys.fes);
    mfem::ConstantCoefficient kappa(130.0);
    sys.a->AddDomainIntegrator(new mfem::DiffusionIntegrator(kappa), attr_mask);
    sys.a->Assemble();

    // ------------------------
    // Load vector: Gaussian
    // ------------------------
    class GaussianSRC : public mfem::Coefficient {
        double cx, cy, cz;
    public:
        GaussianSRC() : cx(0.5), cy(0.5), cz(0.5) {}
        virtual double Eval(mfem::ElementTransformation &T,
                            const mfem::IntegrationPoint &ip)
        {
            mfem::Vector x;
            T.Transform(ip, x);
            double dx = x[0] - cx;
            double dy = x[1] - cy;
            double dz = x[2] - cz;
            return 100.0 * exp(-200.0 * (dx*dx + dy*dy + dz*dz));
        }
    };

    GaussianSRC gsrc;

    sys.b = new mfem::ParLinearForm(sys.fes);
    sys.rhs_int = new mfem::DomainLFIntegrator(gsrc);
    sys.b->AddDomainIntegrator(sys.rhs_int, attr_mask);
    sys.b->Assemble();

    // ------------------------
    // Solution vector
    // ------------------------
    sys.x = new mfem::ParGridFunction(sys.fes);
    *sys.x = 0.0;

    return sys;
}

//--------------------------------------------
// 2. Solve with argument-based solver config
//--------------------------------------------
long SolveSystem(FEMSystem& sys, const std::string& solver_type,
                 bool use_interior_bc, const std::string& prec_type,
                 double rel_tol, int max_iter, int print_level)
{
  // --------------------------------------------------
  // 1) Collect essential true DOFs from boundary attrs
  //    (these are driven by sys.ess_bdr set in Init_FEMSystem)
  // --------------------------------------------------
  mfem::Array<int> bd_tdofs;
  sys.fes->GetEssentialTrueDofs(sys.ess_bdr, bd_tdofs);

  // --------------------------------------------------
  // 2) Optionally add interior DOFs (e.g. x = 0.5 plane)
  // --------------------------------------------------
  mfem::Array<int> active_ess_dofs;

  if (use_interior_bc)
  {
    mfem::Array<int> interior_tdofs;
    // This should fill interior_tdofs with *true dof* indices
    MarkInteriorPlaneDOFs(*sys.pmesh, *sys.fes, interior_tdofs);

    active_ess_dofs = bd_tdofs;
    active_ess_dofs.Append(interior_tdofs);
    active_ess_dofs.Sort();
    active_ess_dofs.Unique();  // IMPORTANT: eliminate duplicates
  }
  else
  {
    active_ess_dofs = bd_tdofs;
  }

  // Sanity check: if we have no essential dofs at all, the system
  // may be singular for a pure Neumann problem with nonzero load.
  if (active_ess_dofs.Size() == 0)
  {
    if (sys.fes->GetParMesh()->GetMyRank() == 0)
      std::cerr << "[SolveSystem] WARNING: No essential DOFs; "
                   "system may be singular for this load.\n";
  }

  // --------------------------------------------------
  // 3) Make sure current solution respects homogeneous BCs
  //    (we assume zero Dirichlet; Init_FEMSystem set *sys.x = 0.0 already)
  //    But we enforce zero again on active essential DOFs just to be safe.
  // --------------------------------------------------
  for (int i = 0; i < active_ess_dofs.Size(); i++)
  {
    (*sys.x)(active_ess_dofs[i]) = 0.0;
  }

  // --------------------------------------------------
  // 4) Form the parallel linear system A X = B with BCs
  // --------------------------------------------------
  mfem::OperatorPtr A;      // owns the operator
  mfem::HypreParVector X;   // solution in true-dof space
  mfem::HypreParVector B;   // RHS in true-dof space

  long residual = -1;

  // This enforces essential BCs encoded in active_ess_dofs and sys.x
  // and modifies A, X, B accordingly.
  sys.a->FormLinearSystem(active_ess_dofs, *sys.x, *sys.b, A, X, B);

  // Extract actual Hypre matrix for Hypre-based preconditioners
  auto *A_hypre = A.As<mfem::HypreParMatrix>();
  MFEM_VERIFY(A_hypre, "FormLinearSystem did not produce a HypreParMatrix.");

  // --------------------------------------------------
  // 5) Build preconditioner
  // --------------------------------------------------
  std::unique_ptr<mfem::Solver> prec;

  if (prec_type == "HypreAMG")
  {
    auto amg = std::make_unique<mfem::HypreBoomerAMG>(*A_hypre);
    prec = std::move(amg);
  }
  else if (prec_type == "Jacobi")
  {
    auto hs = std::make_unique<mfem::HypreSmoother>(*A_hypre);
    hs->SetType(mfem::HypreSmoother::Jacobi);
    prec = std::move(hs);
  }
  else
  {
    if (sys.fes->GetParMesh()->GetMyRank() == 0)
      std::cerr << "[SolveSystem] Unknown preconditioner: "
                << prec_type << std::endl;
    return residual;
  }

  // --------------------------------------------------
  // 6) Choose iterative solver
  // --------------------------------------------------
  std::unique_ptr<mfem::IterativeSolver> solver;
  MPI_Comm comm = sys.fes->GetParMesh()->GetComm();

  if (solver_type == "CG")
  {
    solver = std::make_unique<mfem::CGSolver>(comm);
  }
  else if (solver_type == "MINRES")
  {
    solver = std::make_unique<mfem::MINRESSolver>(comm);
  }
  else if (solver_type == "GMRES")
  {
    solver = std::make_unique<mfem::GMRESSolver>(comm);
  }
  else
  {
    if (sys.fes->GetParMesh()->GetMyRank() == 0)
      std::cerr << "[SolveSystem] Unknown solver: "
                << solver_type << std::endl;
    return residual;
  }

  solver->SetOperator(*A);
  solver->SetPreconditioner(*prec);
  solver->SetRelTol(rel_tol);
  solver->SetAbsTol(0.0);
  solver->SetMaxIter(max_iter);
  solver->SetPrintLevel(print_level);

  // --------------------------------------------------
  // 7) Solve system
  // --------------------------------------------------
  solver->Mult(B, X);

  // If something went wrong (e.g. singular system), norm may be NaN/inf
  double final_norm = solver->GetFinalNorm();
  if (!(final_norm == final_norm) || !std::isfinite(final_norm))  // NaN or inf
  {
    if (sys.fes->GetParMesh()->GetMyRank() == 0)
      std::cerr << "[SolveSystem] WARNING: solver produced non-finite residual "
                << "(likely singular system / bad BCs).\n";
  }

  // --------------------------------------------------
  // 8) Recover FE solution from true-dof solution
  // --------------------------------------------------
  sys.a->RecoverFEMSolution(X, *sys.b, *sys.x);

  if (sys.fes->GetParMesh()->GetMyRank() == 0)
  {
    std::cout << "Solver converged in " << solver->GetNumIterations()
              << " iterations, final residual = " << final_norm
              << std::endl;
  }

  residual = final_norm;
  return residual;
}



} // namespace mfem_support
#endif