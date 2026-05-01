//
// Created by gangwh on 12/1/25.
//
#include <mfem.hpp>
//#include <Omega_h_mesh.hpp>
#ifndef PCMS_MFEM_COUPLING_SUPPORT_H
#define PCMS_MFEM_COUPLING_SUPPORT_H

using namespace mfem;
namespace support
{
/* ------------------------------------------------------------
Writes a ParaView collection for a parallel mesh and a GridFunction solution.
 - pmesh:   parallel mesh (owned elsewhere; must remain valid during Save())
 - sol_gf:  solution GridFunction to visualize (typically defined on pmesh)
 - field_name: name shown in ParaView (e.g., "temperature")
 - base_name:  output collection base (folder/collection name)
 - cycle/time: metadata for ParaView time-series support
------------------------------------------------------------
*/
inline void SaveParaview(ParMesh& pmesh, ParGridFunction& x,
                         const std::string& collection = "thermal_solution",
                         const std::string& field_name = "Temperature",
                         int cycle = 0, double time = 0.0)
{
  ParaViewDataCollection pvdc(collection.c_str(), &pmesh);
  pvdc.SetPrefixPath("paraview");
  pvdc.SetCycle(cycle);
  pvdc.SetTime(time);
  pvdc.SetHighOrderOutput(true);
  pvdc.SetDataFormat(VTKFormat::BINARY);
  pvdc.RegisterField(field_name.c_str(), &x);
  pvdc.Save();
}

struct FEMSystem
{
  mfem::Mesh* mesh = nullptr;
  mfem::ParMesh* pmesh = nullptr;
  mfem::H1_FECollection* fec = nullptr;
  mfem::ParFiniteElementSpace* fes = nullptr;
  mfem::Array<int> ess_tdofs;
  mfem::ParBilinearForm* a = nullptr;
  mfem::ParLinearForm* b = nullptr;
  mfem::ParGridFunction* x = nullptr;
  mfem::DomainLFIntegrator* rhs_int = nullptr;
};

struct ThermalParams
{
  std::array<double, 3> size; // lx, ly, lz
  std::array<int, 3> ne;      // nx, ny, nz
  double q_total;             // total heat [W]
  double kappa;               // thermal conductivity
  double rho;                 // density
  double cp;                  // heat capacity (unused here, steady)
  double h_flux;      // peak heat transfer coefficient for flux BC (used in
                      // FunctionCoefficient)
  double h_conv;      // convection coefficient
  double T_conv;      // ambient temperature for convection BC
  double T_dirichlet; // Dirichlet boundary temperature
};

// -----------------------------
// Utility: tolerance for x comparisons
// -----------------------------
static double DefaultTolX(const ParMesh& pmesh)
{
  double xmin = 1e300, xmax = -1e300;
  for (int i = 0; i < pmesh.GetNV(); i++) {
    const double* v = pmesh.GetVertex(i);
    xmin = std::min(xmin, v[0]);
    xmax = std::max(xmax, v[0]);
  }
  const double Lx = xmax - xmin;
  return std::max(1e-12, 1e-10 * (std::abs(Lx) + 1.0));
}

// -----------------------------
// Mark essential TRUE dofs on boundary line x=a in 2D by boundary elements
// -----------------------------
static void MarkEssTrueDofs_BdrPlaneX_2D(const ParMesh& pmesh,
                                         const ParFiniteElementSpace& fes,
                                         double a, Array<char>& mark,
                                         double tol)
{
  MFEM_VERIFY(pmesh.Dimension() == 2, "This helper is for 2D.");
  Array<int> verts, vdofs;

  for (int be = 0; be < pmesh.GetNBE(); be++) {
    pmesh.GetBdrElementVertices(be, verts);

    bool on_plane = true;
    for (int j = 0; j < verts.Size(); j++) {
      const double* v = pmesh.GetVertex(verts[j]);
      if (std::abs(v[0] - a) > tol) {
        on_plane = false;
        break;
      }
    }
    if (!on_plane) {
      continue;
    }

    fes.GetBdrElementVDofs(be, vdofs);
    for (int k = 0; k < vdofs.Size(); k++) {
      const int tdof = fes.GetLocalTDofNumber(vdofs[k]);
      if (tdof >= 0) {
        mark[tdof] = 1;
      }
    }
  }
}

// -----------------------------
// Mark essential TRUE dofs on internal plane x=a via vertices (mesh-conforming
// only)
// -----------------------------
static void MarkEssTrueDofs_VertPlaneX(const ParMesh& pmesh,
                                       const ParFiniteElementSpace& fes,
                                       double a, Array<char>& mark, double tol)
{
  Array<int> vdofs;
  for (int vi = 0; vi < pmesh.GetNV(); vi++) {
    const double* v = pmesh.GetVertex(vi);
    if (std::abs(v[0] - a) > tol) {
      continue;
    }

    fes.GetVertexVDofs(vi, vdofs);
    for (int k = 0; k < vdofs.Size(); k++) {
      const int tdof = fes.GetLocalTDofNumber(vdofs[k]);
      if (tdof >= 0) {
        mark[tdof] = 1;
      }
    }
  }
}

// -----------------------------
// Convert mark[] -> ess_tdofs list
// -----------------------------
static void MarkToList(const Array<char>& mark, Array<int>& ess_tdofs)
{
  ess_tdofs.SetSize(0);
  for (int i = 0; i < mark.Size(); i++) {
    if (mark[i]) {
      ess_tdofs.Append(i);
    }
  }
}

// -----------------------------
// Initialize Dirichlet values in sys.x (TRUE dof vector) for elimination.
// For order=1 nodal H1, vertex dofs are enough.
// -----------------------------
static void InitializeDirichletValues_Order1(const ParMesh& pmesh,
                                             const ParFiniteElementSpace& fes,
                                             const Array<char>& mark,
                                             ParGridFunction& x, double tol,
                                             double T_left, double T_right,
                                             double T_mid)
{
  Vector xt(fes.GetTrueVSize());
  x.GetTrueDofs(xt);

  Array<int> vdofs;

  for (int vi = 0; vi < pmesh.GetNV(); vi++) {
    const double* v = pmesh.GetVertex(vi);
    const double xv = v[0];

    const bool is_left = (std::abs(xv - 0.0) <= tol);
    const bool is_right = (std::abs(xv - 1.0) <= tol);
    const bool is_mid = (std::abs(xv - 0.5) <= tol);

    if (!(is_left || is_right || is_mid)) {
      continue;
    }

    double val = 0.0;
    if (is_left)
      val = T_left;
    else if (is_right)
      val = T_right;
    else
      val = T_mid;

    fes.GetVertexVDofs(vi, vdofs);
    for (int k = 0; k < vdofs.Size(); k++) {
      const int tdof = fes.GetLocalTDofNumber(vdofs[k]);
      if (tdof < 0) {
        continue;
      }
      if (!mark[tdof]) {
        continue;
      }
      xt[tdof] = val;
    }
  }

  x.SetFromTrueDofs(xt);
}

// -----------------------------
// Build FE system: -div(k grad T) = 0
// -----------------------------
inline FEMSystem  Init_FEMSystem(ParMesh* pmesh, int order, double kappa_val)
{
  FEMSystem sys;
  sys.pmesh = pmesh;
  const int dim = pmesh->Dimension();

  sys.fec = new H1_FECollection(order, dim);
  sys.fes = new ParFiniteElementSpace(pmesh, sys.fec);

  ConstantCoefficient zero(0.0);
  sys.b = new ParLinearForm(sys.fes);
  sys.b->AddDomainIntegrator(new DomainLFIntegrator(zero));
  sys.b->Assemble();

  ConstantCoefficient kappa(kappa_val);
  sys.a = new ParBilinearForm(sys.fes);
  sys.a->AddDomainIntegrator(new DiffusionIntegrator(kappa));
  sys.a->Assemble();
  sys.a->Finalize();

  sys.x = new ParGridFunction(sys.fes);
  *sys.x = 0.0;

  return sys;
}
inline void  DestroyFEMSystem(FEMSystem& sys)
{
  delete sys.x;
  delete sys.a;
  delete sys.b;
  delete sys.fes;
  delete sys.fec;

  sys.x = nullptr;
  sys.a = nullptr;
  sys.b = nullptr;
  sys.fes = nullptr;
  sys.fec = nullptr;
  sys.pmesh = nullptr;
}
// -----------------------------
// SolveSystem (low verbosity)
// -----------------------------
inline long SolveSystem(FEMSystem& sys, const std::string& solver_type,
                 const std::string& prec_type, double rel_tol, int max_iter)
{
  OperatorPtr A;
  HypreParVector X, B;

  sys.a->FormLinearSystem(sys.ess_tdofs, *sys.x, *sys.b, A, X, B);
  auto* A_hypre = A.As<HypreParMatrix>();
  MFEM_VERIFY(A_hypre, "FormLinearSystem did not produce HypreParMatrix.");

  std::unique_ptr<Solver> prec;
  if (prec_type == "HypreAMG") {
    auto amg = std::make_unique<HypreBoomerAMG>(*A_hypre);
    amg->SetPrintLevel(0);
    prec = std::move(amg);
  } else if (prec_type == "Jacobi") {
    auto sm = std::make_unique<HypreSmoother>(*A_hypre);
    sm->SetType(HypreSmoother::Jacobi);
    prec = std::move(sm);
  } else {
    MFEM_ABORT("Unknown preconditioner.");
  }

  std::unique_ptr<IterativeSolver> solver;
  MPI_Comm comm = sys.fes->GetParMesh()->GetComm();

  if (solver_type == "CG")
    solver = std::make_unique<CGSolver>(comm);
  else if (solver_type == "MINRES")
    solver = std::make_unique<MINRESSolver>(comm);
  else if (solver_type == "GMRES")
    solver = std::make_unique<GMRESSolver>(comm);
  else
    MFEM_ABORT("Unknown solver.");

  solver->SetOperator(*A_hypre);
  solver->SetPreconditioner(*prec);
  solver->SetRelTol(rel_tol);
  solver->SetAbsTol(0.0);
  solver->SetMaxIter(max_iter);
  solver->SetPrintLevel(0);

  solver->Mult(B, X);
  sys.a->RecoverFEMSolution(X, *sys.b, *sys.x);

  return solver->GetFinalNorm();
}

// -----------------------------
// Diagnostics: count constrained vertices, and report Tmin/Tmax on x=0.5 line
// -----------------------------
static void ReportLineStats_Order1(const ParMesh& pmesh,
                                   const ParGridFunction& x, double xline,
                                   double tol, const char* label)
{
  double Tmin = 1e300, Tmax = -1e300;
  int count = 0;

  for (int vi = 0; vi < pmesh.GetNV(); vi++) {
    const double* v = pmesh.GetVertex(vi);
    if (std::abs(v[0] - xline) <= tol) {
      const double Tv = x(vi); // order=1 nodal H1 -> vertex value access
      Tmin = std::min(Tmin, Tv);
      Tmax = std::max(Tmax, Tv);
      count++;
    }
  }

  std::cout << label << ": vertices=" << count << "  Tmin=" << Tmin
            << "  Tmax=" << Tmax << "\n";
}
} // namespace support
#endif // PCMS_MFEM_COUPLING_SUPPORT_H
