
// Notes:
// - Mesh: unit square [0,1]x[0,1] with 10x10 QUADs
// - Thermal BC: Dirichlet T = T_left on x=0; insulated elsewhere
// - Elastic BC: clamped (u=0) on x=0; free elsewhere
// - Coupling: elasticity RHS gets body force f = alpha*(lambda*dim +
// 2*mu)*grad(T)
// - Loosely-coupled iteration stops when relative L2 changes in T and u fall
// below tol

#include "mfem.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace mfem;
using namespace std;

// Mark boundary attributes by geometric side for a unit square mesh.
// 1: left (x=0), 2: right (x=1), 3: bottom (y=0), 4: top (y=1)
static void MarkUnitSquareSides(Mesh& mesh)
{
  const int dim = mesh.Dimension();
  MFEM_VERIFY(dim == 2, "This helper expects a 2D unit square mesh.");
  Vector xc(dim);
  for (int be = 0; be < mesh.GetNBE(); ++be) {
    Element* el = mesh.GetBdrElement(be);
    Array<int> v;
    el->GetVertices(v);
    xc = 0.0;
    for (int j = 0; j < v.Size(); ++j) {
      const double* X = mesh.GetVertex(v[j]);
      xc[0] += X[0];
      xc[1] += X[1];
    }
    xc /= v.Size();
    int attr = 0;
    if (fabs(xc[0] - 0.0) < 1e-8)
      attr = 1; // left
    else if (fabs(xc[0] - 1.0) < 1e-8)
      attr = 2; // right
    else if (fabs(xc[1] - 0.0) < 1e-8)
      attr = 3; // bottom
    else if (fabs(xc[1] - 1.0) < 1e-8)
      attr = 4; // top
    else
      attr = 99; // unknown (should not happen on unit square)
    mesh.SetBdrAttribute(be, attr);
  }
}

// Solve: -div(kappa grad T) = 0 with Dirichlet BCs encoded in ess_bdr_T
void SolveThermal(FiniteElementSpace& T_fes, GridFunction& T_gf,
                  const Array<int>& ess_bdr_T, double kappa, double T_left,
                  double T_ref)
{
  ConstantCoefficient kappa_c(kappa);
  ConstantCoefficient zero(0.0);
  ConstantCoefficient T_left_c(T_left);

  BilinearForm aT(&T_fes);
  aT.AddDomainIntegrator(new DiffusionIntegrator(kappa_c));

  LinearForm bT(&T_fes);
  bT.AddDomainIntegrator(new DomainLFIntegrator(zero)); // zero source

  aT.Assemble();
  bT.Assemble();

  // Apply Dirichlet BC on selected boundaries
  Array<int> ess_tdof_list_T;
  T_fes.GetEssentialTrueDofs(ess_bdr_T, ess_tdof_list_T);
  // Project boundary values
  T_gf.ProjectBdrCoefficient(T_left_c, ess_bdr_T);

  OperatorPtr AopT;
  Vector B_T, X_T;
  aT.FormLinearSystem(ess_tdof_list_T, T_gf, bT, AopT, X_T, B_T);

  CGSolver cg;
  cg.SetRelTol(1e-12);
  cg.SetAbsTol(0.0);
  cg.SetMaxIter(5000);
  cg.SetPrintLevel(0);
  GSSmoother M((SparseMatrix&)(*AopT));
  cg.SetOperator(*AopT);
  cg.SetPreconditioner(M);
  cg.Mult(B_T, X_T);

  aT.RecoverFEMSolution(X_T, bT, T_gf);
}

// Solve linear elasticity with thermal-expansion body force f(T) = Kth *
// grad(T) Ess. BCs (clamped) encoded in ess_bdr_U
void SolveElastic(FiniteElementSpace& U_fes, GridFunction& U_gf,
                  const Array<int>& ess_bdr_U, const GridFunction& T_gf,
                  double lambda, double mu, double alpha)
{
  const int dim = U_fes.GetMesh()->Dimension();
  const double Kth = alpha * (lambda * dim + 2.0 * mu);

  ConstantCoefficient lambda_c(lambda), mu_c(mu);

  // --- Elastic stiffness ---
  BilinearForm aU(&U_fes);
  aU.AddDomainIntegrator(new ElasticityIntegrator(lambda_c, mu_c));

  // --- Thermal body force ---
  GradientGridFunctionCoefficient gradT(&T_gf);
  auto fth = gradT * Kth; // this returns a temporary scaled VectorCoefficient
  LinearForm bU(&U_fes);
  bU.AddDomainIntegrator(new VectorDomainLFIntegrator(fth));

  // --- Assemble system ---
  aU.Assemble();
  bU.Assemble();

  Array<int> ess_tdof_list_U;
  U_fes.GetEssentialTrueDofs(ess_bdr_U, ess_tdof_list_U);

  OperatorPtr AopU;
  Vector B_U, X_U;
  aU.FormLinearSystem(ess_tdof_list_U, U_gf, bU, AopU, X_U, B_U);

  CGSolver cg;
  cg.SetRelTol(1e-12);
  cg.SetAbsTol(0.0);
  cg.SetMaxIter(10000);
  cg.SetPrintLevel(0);

  GSSmoother M((SparseMatrix&)(*AopU));
  cg.SetOperator(*AopU);
  cg.SetPreconditioner(M);
  cg.Mult(B_U, X_U);

  aU.RecoverFEMSolution(X_U, bU, U_gf);
}

int main(int argc, char* argv[])
{
  // --- Parameters ---
  const int nx = 10, ny = 10;
  const int order = 1;
  const double kappa = 45.0; // W/m-K
  const double E = 200e9;    // Pa
  const double nu = 0.30;
  const double alpha = 1.0e-5;         // 1/K
  const double T_ref = 293.15;         // K
  const double T_left = T_ref + 100.0; // K

  const double lambda = (E * nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
  const double mu = E / (2.0 * (1.0 + nu));

  const int max_iter = 50;
  const double tol_T = 1e-10;
  const double tol_U = 1e-10;
  const double omega_T = 1.0; // under-relaxation (if desired)
  const double omega_U = 1.0;

  // --- Mesh & spaces ---
  Mesh mesh(nx, ny, Element::QUADRILATERAL, true, 1.0, 1.0);
  MarkUnitSquareSides(mesh);
  mesh.EnsureNodes();

  const int dim = mesh.Dimension();
  H1_FECollection h1_fec(order, dim);
  FiniteElementSpace T_fes(&mesh, &h1_fec);      // scalar T
  FiniteElementSpace U_fes(&mesh, &h1_fec, dim); // vector u

  // --- Essential boundary markers ---
  Array<int> ess_bdr_T(mesh.bdr_attributes.Max());
  ess_bdr_T = 0; // Dirichlet for T
  Array<int> ess_bdr_U(mesh.bdr_attributes.Max());
  ess_bdr_U = 0; // clamp for u
  if (mesh.bdr_attributes.Size() > 0) {
    ess_bdr_T[0] = 1; // attribute 1: left
    ess_bdr_U[0] = 1; // attribute 1: left
  }

  // --- Unknowns ---
  GridFunction T_gf(&T_fes);
  T_gf = T_ref;
  GridFunction U_gf(&U_fes);
  U_gf = 0.0;

  GridFunction T_old(&T_fes);
  T_old = T_gf;
  GridFunction U_old(&U_fes);
  U_old = U_gf;

  // --- Coupling loop ---
  cout << fixed << setprecision(3);
  for (int it = 0; it < max_iter; ++it) {
    // Store previous
    T_old = T_gf;
    U_old = U_gf;

    // Thermal solve (independent in this model)
    SolveThermal(T_fes, T_gf, ess_bdr_T, kappa, T_left, T_ref);

    // Under-relax T if needed
    if (omega_T != 1.0) {
      T_gf *= omega_T;
      T_gf.Add(1.0 - omega_T, T_old);
    }

    // Elastic solve (depends on T)
    SolveElastic(U_fes, U_gf, ess_bdr_U, T_gf, lambda, mu, alpha);

    // Under-relax U if needed
    if (omega_U != 1.0) {
      U_gf *= omega_U;
      U_gf.Add(1.0 - omega_U, U_old);
    }

    // Convergence checks (relative L2 changes)
    double dT = T_gf.DistanceTo(T_old); // ||T - Told||
    double nT = T_gf.Norml2();          // ||T||
    double relT = dT / (nT + 1e-30);

    double dU = U_gf.DistanceTo(U_old); // ||u - uold||
    double nU = U_gf.Norml2();          // ||u||
    double relU = dU / (nU + 1e-30);

    cout << "Iter " << setw(3) << it << "  relT=" << scientific << relT
         << "  relU=" << scientific << relU << defaultfloat << endl;

    if (relT < tol_T && relU < tol_U) {
      cout << "Converged at iteration " << it << endl;
      break;
    }

    if (it == max_iter - 1) {
      cout << "Reached max_iter without meeting tolerances." << endl;
    }
  }

  // --- Output ---
  {
    ofstream ot("T.gf");
    T_gf.Save(ot);
    ofstream ou("U.gf");
    U_gf.Save(ou);
  }

  VisItDataCollection visit_dc("thermo_elastic_serial_iter", &mesh);
  visit_dc.RegisterField("temperature", &T_gf);
  visit_dc.RegisterField("displacement", &U_gf);
  visit_dc.SetCycle(0);
  visit_dc.SetTime(0.0);
  visit_dc.Save();

  cout << "Done. Wrote VisIt files in ./thermo_elastic_serial_iter/ and "
          "GridFunctions T.gf, U.gf\n";
  return 0;
}
