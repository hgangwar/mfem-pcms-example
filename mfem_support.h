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

//--------------------------------------------
// 1. Initialization (Mesh + FE + BilinearForm)
//--------------------------------------------
FEMSystem Init_FEMSystem(MPI_Comm comm, const std::string& mesh_file, int order,
                         const char client)
{
  FEMSystem sys;

  sys.mesh = new mfem::Mesh(mesh_file.c_str(), 1, 1, true);
  sys.pmesh = new mfem::ParMesh(comm, *sys.mesh);

  int dim = sys.pmesh->Dimension();

  // ---------------------------
  // Assign element attributes
  // ---------------------------
  AssignAttributesByX(*sys.pmesh);

  // ---------------------------
  // Build attr_mask based on client type
  // (element attributes 1,2,3 from AssignAttributesByX)
  // ---------------------------
  const int max_attr = sys.pmesh->attributes.Max(); // should be 3
  mfem::Array<int> attr_mask(max_attr);
  attr_mask = 0;

  switch (client) {
    case 'A':
      // client A: region 1 + overlap (1,2)
      attr_mask[0] = 1; // attribute 1
      attr_mask[1] = 1; // attribute 2
      break;
    case 'B':
      // client B: overlap + region 3 (2,3)
      attr_mask[1] = 1; // attribute 2
      attr_mask[2] = 1; // attribute 3
      break;
    default: throw std::invalid_argument("Unknown client type");
  }

  // ---------------------------
  // Boundary attribute mask for external Dirichlet BCs
  // ---------------------------
  sys.ess_bdr.SetSize(sys.pmesh->bdr_attributes.Max());
  sys.ess_bdr = 0;
  for (int k = 0; k < sys.ess_bdr.Size(); k++)
    sys.ess_bdr[k] = 1;

  // FE space
  sys.fec = new mfem::H1_FECollection(order, dim);
  sys.fes = new mfem::ParFiniteElementSpace(sys.pmesh, sys.fec);

  // Bilinear form (domain restricted by attr_mask)
  sys.a = new mfem::ParBilinearForm(sys.fes);
  mfem::ConstantCoefficient kappa(130.0);
  sys.a->AddDomainIntegrator(new mfem::DiffusionIntegrator(kappa), attr_mask);
  sys.a->Assemble();
  sys.a->Finalize();

  // Linear form (same attr_mask)
  mfem::ConstantCoefficient f(100.0);
  sys.b = new mfem::ParLinearForm(sys.fes);
  sys.rhs_int = new mfem::DomainLFIntegrator(f);
  sys.b->AddDomainIntegrator(sys.rhs_int, attr_mask);
  sys.b->Assemble();

  // Solution
  sys.x = new mfem::ParGridFunction(sys.fes);
  *sys.x = 0.0;

  return sys;
}

//--------------------------------------------
// 2. Solve with argument-based solver config
//--------------------------------------------
long SolveSystem(FEMSystem& sys, const std::string& solver_type,
                 bool use_interior_bc, const std::string& prec_type,
                 double rel_tol = 1e-8, int max_iter = 500, int print_level = 5)
{
  // Boundary essential DOFs from boundary attributes
  mfem::Array<int> bd_tdofs;
  sys.fes->GetEssentialTrueDofs(sys.ess_bdr, bd_tdofs);

  // Interior DOFs on x = 0.5
  mfem::Array<int> active_ess_dofs;

  if (use_interior_bc) {
    mfem::Array<int> interior_tdofs;
    MarkInteriorPlaneDOFs(*sys.pmesh, *sys.fes, interior_tdofs);
    // union: boundary + interior
    active_ess_dofs = bd_tdofs;
    active_ess_dofs.Append(interior_tdofs);
    active_ess_dofs.Sort(); // optional: sort + dedup
  } else {
    active_ess_dofs = bd_tdofs;
  }

  // 3) Form the parallel linear system A X = B
  mfem::OperatorPtr A; // Owned smart pointer
  mfem::HypreParVector X, B;
  long residual = -1;
  sys.a->FormLinearSystem(active_ess_dofs, *sys.x, *sys.b, A, X, B);

  // Extract actual Hypre matrix pointer for Hypre-based preconditioners
  auto* A_hypre = A.As<mfem::HypreParMatrix>();
  MFEM_VERIFY(A_hypre, "FormLinearSystem did not produce a HypreParMatrix.");

  // --------------------------------------------------
  // Select preconditioner (Hypre-only)
  // --------------------------------------------------
  std::unique_ptr<mfem::Solver> prec;

  if (prec_type == "HypreAMG") {
    auto amg = std::make_unique<mfem::HypreBoomerAMG>(*A_hypre);
    prec = std::move(amg);
  } else if (prec_type == "Jacobi") {
    auto hs = std::make_unique<mfem::HypreSmoother>(*A_hypre);
    hs->SetType(mfem::HypreSmoother::Jacobi);
    prec = std::move(hs);
  } else {
    if (sys.fes->GetParMesh()->GetMyRank() == 0)
      std::cerr << "Unknown preconditioner: " << prec_type << std::endl;
    return residual;
  }

  // --------------------------------------------------
  // Choose parallel iterative solver
  // --------------------------------------------------
  std::unique_ptr<mfem::IterativeSolver> solver;
  MPI_Comm comm = sys.fes->GetParMesh()->GetComm();

  if (solver_type == "CG")
    solver = std::make_unique<mfem::CGSolver>(comm);
  else if (solver_type == "MINRES")
    solver = std::make_unique<mfem::MINRESSolver>(comm);
  else if (solver_type == "GMRES")
    solver = std::make_unique<mfem::GMRESSolver>(comm);
  else {
    if (sys.fes->GetParMesh()->GetMyRank() == 0)
      std::cerr << "Unknown solver: " << solver_type << std::endl;
    return residual;
  }
  solver->SetOperator(*A);
  solver->SetPreconditioner(*prec);
  solver->SetRelTol(rel_tol);
  solver->SetAbsTol(0.0);
  solver->SetMaxIter(max_iter);
  solver->SetPrintLevel(print_level);

  // --------------------------------------------------
  // Solve system
  // --------------------------------------------------
  solver->Mult(B, X);

  // --------------------------------------------------
  // Recover finite element solution
  // --------------------------------------------------
  sys.a->RecoverFEMSolution(X, *sys.b, *sys.x);

  if (sys.fes->GetParMesh()->GetMyRank() == 0) {
    std::cout << "Solver converged in " << solver->GetNumIterations()
              << " iterations, final residual = " << solver->GetFinalNorm()
              << std::endl;
  }
  residual = solver->GetFinalNorm();
  return residual;
}

} // namespace mfem_support
#endif