//
// Created by gangwh on 12/1/25.
//
#include<mfem.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_mesh.hpp>
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
inline void SaveParaviewSolution(mfem::Mesh &mesh,
                                mfem::GridFunction &sol_gf,
                                const std::string &field_name = "solution",
                                const std::string &base_name  = "paraview_output",
                                int cycle = 0,
                                double time = 0.0,
                                bool high_order = true)
{
  mfem::ParaViewDataCollection pvd(base_name.c_str(), &mesh);
  pvd.RegisterField(field_name.c_str(), &sol_gf);
  pvd.SetCycle(cycle);
  pvd.SetTime(time);
  pvd.SetHighOrderOutput(high_order);

  pvd.Save();
}

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

  std::array<double,3> size;  // domain size: Lx, Ly, Lz
  char client = 'M';          // 'A', 'B', or 'M'
};
struct ThermalParams
{
  std::array<double,3> size;      // lx, ly, lz
  std::array<int,3>    ne;        // nx, ny, nz
  double q_total;                 // total heat [W]
  double kappa;                   // thermal conductivity
  double rho;                     // density
  double cp;                      // heat capacity (unused here, steady)
  double h_flux;                  // peak heat transfer coefficient for flux BC (used in FunctionCoefficient)
  double h_conv;                  // convection coefficient
  double T_conv;                  // ambient temperature for convection BC
  double T_dirichlet;             // Dirichlet boundary temperature
};

// ------------------------------------------------------------
// Utilities: attribute assignment by x
// ------------------------------------------------------------
void AssignAttributesByX(ParMesh &pmesh, double Lx)
{
  const int ne = pmesh.GetNE();
  for (int e = 0; e < ne; e++)
  {
    Element *el = pmesh.GetElement(e);
    const int *verts = el->GetVertices();
    int nv = el->GetNVertices();

    double xc = 0.0;
    for (int j = 0; j < nv; j++)
    {
      const double *coord = pmesh.GetVertex(verts[j]);
      xc += coord[0];
    }
    xc /= nv;

    int attr;
    if (xc < 0.4 * Lx)       attr = 1;
    else if (xc <= 0.6 * Lx) attr = 2;
    else                     attr = 3;

    el->SetAttribute(attr);
  }
}

// ------------------------------------------------------------
// Initialize x with external Dirichlet values (300 at x=0, 350 at x=Lx)
// ------------------------------------------------------------
void ApplyExternalDirichlet(const ParMesh &pmesh,
                            ParFiniteElementSpace &fes,
                            ParGridFunction &x,
                            double Lx)
{
  const double eps = 1e-6;
  Array<int> vdofs;

  for (int v = 0; v < pmesh.GetNV(); v++)
  {
    const double *coord = pmesh.GetVertex(v);
    double X = coord[0];

    bool on_left  = (X < eps);
    bool on_right = (X > Lx - eps);

    if (!on_left && !on_right) continue;

    fes.GetVertexDofs(v, vdofs);
    for (int k = 0; k < vdofs.Size(); k++)
    {
      int tdof = vdofs[k];
      if (tdof < 0) tdof = -1 - tdof;

      double T = on_left ? 300.0 : 350.0;
      x(tdof) = T;
    }
  }
}

// ------------------------------------------------------------
// Mark external Dirichlet DOFs for essential list
// ------------------------------------------------------------
void MarkExternalDirichletDOFs(const ParMesh &pmesh,
                               const ParFiniteElementSpace &fes,
                               double Lx,
                               Array<int> &tdofs_out)
{
  const double eps = 1e-6;
  Array<int> vdofs;
  Array<int> marker(fes.GetTrueVSize());
  marker = 0;

  for (int v = 0; v < pmesh.GetNV(); v++)
  {
    const double *coord = pmesh.GetVertex(v);
    double X = coord[0];

    bool on_left  = (X < eps);
    bool on_right = (X > Lx - eps);

    if (!on_left && !on_right) continue;

    fes.GetVertexDofs(v, vdofs);
    for (int k = 0; k < vdofs.Size(); k++)
    {
      int td = vdofs[k];
      if (td < 0) td = -1 - td;
      marker[td] = 1;
    }
  }

  tdofs_out.DeleteAll();
  for (int i = 0; i < marker.Size(); i++)
    if (marker[i]) tdofs_out.Append(i);
}

// ------------------------------------------------------------
// Mark DOFs on internal interface x = 0.5 * Lx
// ------------------------------------------------------------
void MarkInternalInterfaceDOFs(const ParMesh &pmesh,
                               const ParFiniteElementSpace &fes,
                               double Lx,
                               Array<int> &tdofs_out)
{
  const double eps = 1e-6;
  double Xmid = 0.5 * Lx;

  Array<int> vdofs;
  Array<int> marker(fes.GetTrueVSize());
  marker = 0;

  for (int v = 0; v < pmesh.GetNV(); v++)
  {
    const double *coord = pmesh.GetVertex(v);
    double X = coord[0];

    bool on_mid = (fabs(X - Xmid) < eps);
    if (!on_mid) continue;

    fes.GetVertexDofs(v, vdofs);
    for (int k = 0; k < vdofs.Size(); k++)
    {
      int td = vdofs[k];
      if (td < 0) td = -1 - td;
      marker[td] = 1;
    }
  }

  tdofs_out.DeleteAll();
  for (int i = 0; i < marker.Size(); i++)
    if (marker[i]) tdofs_out.Append(i);
}

// ------------------------------------------------------------
// Clamp inactive region temperatures in x to 273 K (for clients A,B)
//   A: inactive region 3 → x > 0.6 * Lx
//   B: inactive region 1 → x < 0.4 * Lx
//   M: no inactive region → no clamping
// ------------------------------------------------------------
void ClampInactiveRegions(FEMSystem &sys)
{
  const double Lx  = sys.size[0];
  const double eps = 1e-6;
  const double Tinactive = 273.0;

  if (sys.client == 'M') return; // all regions active

  Array<int> vdofs;

  for (int v = 0; v < sys.pmesh->GetNV(); v++)
  {
    const double *coord = sys.pmesh->GetVertex(v);
    double X = coord[0];

    bool in_inactive = false;

    if (sys.client == 'A') // active 1,2 → inactive region is 3 (x > 0.6 Lx)
    {
      in_inactive = (X > 0.6 * Lx + eps);
    }
    else if (sys.client == 'B') // active 2,3 → inactive region is 1 (x < 0.4 Lx)
    {
      in_inactive = (X < 0.4 * Lx - eps);
    }

    if (!in_inactive) continue;

    sys.fes->GetVertexDofs(v, vdofs);
    for (int k = 0; k < vdofs.Size(); k++)
    {
      int tdof = vdofs[k];
      if (tdof < 0) tdof = -1 - tdof;
      (*sys.x)(tdof) = Tinactive;
    }
  }
}



// ------------------------------------------------------------
// Init_FEMSystem with low-k inactive regions & heat in region 2
// ------------------------------------------------------------
FEMSystem Init_FEMSystem(MPI_Comm comm,
                         const ThermalParams &params,
                         int order,
                         char client)
{
  FEMSystem sys;
  sys.size   = params.size;
  sys.client = client;

  Mesh serial =
    Mesh::MakeCartesian3D(params.ne[0], params.ne[1], params.ne[2],
                          Element::HEXAHEDRON,
                          params.size[0], params.size[1], params.size[2],
                          true);

  sys.mesh  = new Mesh(serial);
  sys.pmesh = new ParMesh(comm, *sys.mesh);

  // Set attributes 1,2,3 by x
  AssignAttributesByX(*sys.pmesh, params.size[0]);

  int dim = sys.pmesh->Dimension();

  // Active attributes by client
  const int max_attr = sys.pmesh->attributes.Max();
  Array<int> attr_active(max_attr);
  attr_active = 0;

  switch (client)
  {
    case 'A': attr_active[0] = 1; attr_active[1] = 1; break;  // 1,2 active
    case 'B': attr_active[1] = 1; attr_active[2] = 1; break;  // 2,3 active
    case 'M': for (int i = 0; i < max_attr; i++) attr_active[i] = 1; break;
    default:  MFEM_ABORT("Unknown client");
  }

  // Attributes where heat is applied: only region 2
  Array<int> attr_heat(max_attr);
  attr_heat = 0;
  if (max_attr >= 2) attr_heat[1] = 1; // attribute=2 → index 1

  sys.fec = new H1_FECollection(order, dim);
  sys.fes = new ParFiniteElementSpace(sys.pmesh, sys.fec);

  // Volume heat density (uniform within region 2)
  double volume    = params.size[0]*params.size[1]*params.size[2];
  double q_vol_val = params.q_total / volume;
  ConstantCoefficient q_vol(q_vol_val);

  sys.b = new ParLinearForm(sys.fes);
  sys.b->AddDomainIntegrator(new DomainLFIntegrator(q_vol), attr_heat);
  sys.b->Assemble();

  // Piecewise kappa: active → kappa, inactive → k_inactive
  const double k_inactive = 1e-12;
  Vector kvals(max_attr);
  for (int a = 0; a < max_attr; a++)
  {
    if (attr_active[a]) kvals[a] = params.kappa;
    else                kvals[a] = k_inactive;
  }
  PWConstCoefficient kappa_piece(kvals);

  sys.a = new ParBilinearForm(sys.fes);
  sys.a->AddDomainIntegrator(new DiffusionIntegrator(kappa_piece));
  sys.a->Assemble();
  sys.a->Finalize();

  sys.x = new ParGridFunction(sys.fes);
  *sys.x = 0.0;

  // Initialize external BCs in x: 300 at x=0, 350 at x=Lx
  ApplyExternalDirichlet(*sys.pmesh, *sys.fes, *sys.x, params.size[0]);

  int nbdr = sys.pmesh->bdr_attributes.Max();
  sys.ess_bdr.SetSize(nbdr);
  sys.ess_bdr = 0;  // we do geometry-based BCs

  return sys;
}
// ------------------------------------------------------------
// SolveSystem with optional internal BC at x=0.5*Lx.
//   use_internal_bc = false → only external Dirichlet.
//   use_internal_bc = true  → also impose T(x=0.5) from sys.x.
// ------------------------------------------------------------
long SolveSystem(FEMSystem &sys,
                 const std::string &solver_type,
                 bool use_internal_bc,
                 const std::string &prec_type,
                 double rel_tol,
                 int max_iter,
                 int print_level)
{
  Array<int> bd_tdofs;
  MarkExternalDirichletDOFs(*sys.pmesh, *sys.fes, sys.size[0], bd_tdofs);

  Array<int> active_ess = bd_tdofs;

  if (use_internal_bc)
  {
    Array<int> mid_tdofs;
    MarkInternalInterfaceDOFs(*sys.pmesh, *sys.fes, sys.size[0], mid_tdofs);
    active_ess.Append(mid_tdofs);
    active_ess.Sort();
    active_ess.Unique();
    // sys.x already holds the "guess" on x=0.5 from previous solve/coupler.
  }

  OperatorPtr A;
  HypreParVector X, B;

  sys.a->FormLinearSystem(active_ess, *sys.x, *sys.b, A, X, B);

  auto *A_hypre = A.As<HypreParMatrix>();
  MFEM_VERIFY(A_hypre, "FormLinearSystem did not produce a HypreParMatrix.");

  // Preconditioner
  std::unique_ptr<Solver> prec;
  if (prec_type == "HypreAMG")
  {
    auto amg = std::make_unique<HypreBoomerAMG>(*A_hypre);
    prec = std::move(amg);
  }
  else if (prec_type == "Jacobi")
  {
    auto hs = std::make_unique<HypreSmoother>(*A_hypre);
    hs->SetType(HypreSmoother::Jacobi);
    prec = std::move(hs);
  }
  else
  {
    if (sys.fes->GetParMesh()->GetMyRank() == 0)
      std::cerr << "Unknown preconditioner: " << prec_type << std::endl;
    return -1;
  }

  // Solver
  std::unique_ptr<IterativeSolver> solver;
  MPI_Comm comm = sys.fes->GetParMesh()->GetComm();

  if (solver_type == "CG")
    solver = std::make_unique<CGSolver>(comm);
  else if (solver_type == "MINRES")
    solver = std::make_unique<MINRESSolver>(comm);
  else if (solver_type == "GMRES")
    solver = std::make_unique<GMRESSolver>(comm);
  else
  {
    if (sys.fes->GetParMesh()->GetMyRank() == 0)
      std::cerr << "Unknown solver: " << solver_type << std::endl;
    return -1;
  }

  solver->SetOperator(*A);
  solver->SetPreconditioner(*prec);
  solver->SetRelTol(rel_tol);
  solver->SetAbsTol(0.0);
  solver->SetMaxIter(max_iter);
  solver->SetPrintLevel(print_level);

  solver->Mult(B, X);

  sys.a->RecoverFEMSolution(X, *sys.b, *sys.x);

  // Clamp inactive regions for output (A/B clients only)
  ClampInactiveRegions(sys);

  long residual = solver->GetFinalNorm();

  if (sys.fes->GetParMesh()->GetMyRank() == 0)
  {
    std::cout << "Solver converged in " << solver->GetNumIterations()
              << " iterations, final residual = " << residual << std::endl;
  }

  return residual;
}
}
#endif // PCMS_MFEM_COUPLING_SUPPORT_H
