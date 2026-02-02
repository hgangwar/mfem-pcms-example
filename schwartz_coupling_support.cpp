//
// Created by gangwh on 2/2/26.
//

#include "mfem.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <map>
#include <iomanip>

using namespace mfem;

// -----------------------------
// FEMSystem
// -----------------------------
struct FEMSystem
{
  ParMesh *pmesh = nullptr;
  H1_FECollection *fec = nullptr;
  ParFiniteElementSpace *fes = nullptr;

  ParBilinearForm *a = nullptr;
  ParLinearForm *b = nullptr;
  ParGridFunction *x = nullptr;

  Array<int> ess_tdofs; // essential TRUE dofs for elimination
};

// -----------------------------
// Build FE system: -div(k grad T)=0
// -----------------------------
static FEMSystem InitThermalSystem(ParMesh *pmesh, int order, double kappa_val)
{
  FEMSystem sys;
  sys.pmesh = pmesh;

  const int dim = pmesh->Dimension();
  sys.fec = new H1_FECollection(order, dim);
  sys.fes = new ParFiniteElementSpace(pmesh, sys.fec);

  // RHS = 0
  ConstantCoefficient zero(0.0);
  sys.b = new ParLinearForm(sys.fes);
  sys.b->AddDomainIntegrator(new DomainLFIntegrator(zero));
  sys.b->Assemble();

  // diffusion
  ConstantCoefficient kappa(kappa_val);
  sys.a = new ParBilinearForm(sys.fes);
  sys.a->AddDomainIntegrator(new DiffusionIntegrator(kappa));
  sys.a->Assemble();
  sys.a->Finalize();

  sys.x = new ParGridFunction(sys.fes);
  *sys.x = 0.0;

  return sys;
}
// -----------------------------
// SolveSystem
// -----------------------------
static double SolveSystem(FEMSystem &sys,
                          const std::string &solver_type,
                          const std::string &prec_type,
                          double rel_tol,
                          int max_iter)
{
  OperatorPtr A;
  HypreParVector X, B;

  sys.a->FormLinearSystem(sys.ess_tdofs, *sys.x, *sys.b, A, X, B);

  auto *A_hypre = A.As<HypreParMatrix>();
  MFEM_VERIFY(A_hypre, "FormLinearSystem did not produce HypreParMatrix.");

  std::unique_ptr<Solver> prec;
  if (prec_type == "HypreAMG")
  {
    auto amg = std::make_unique<HypreBoomerAMG>(*A_hypre);
    amg->SetPrintLevel(0);
    prec = std::move(amg);
  }
  else if (prec_type == "Jacobi")
  {
    auto sm = std::make_unique<HypreSmoother>(*A_hypre);
    sm->SetType(HypreSmoother::Jacobi);
    prec = std::move(sm);
  }
  else
  {
    MFEM_ABORT("Unknown preconditioner.");
  }

  std::unique_ptr<IterativeSolver> solver;
  MPI_Comm comm = sys.fes->GetParMesh()->GetComm();

  if (solver_type == "CG")          solver = std::make_unique<CGSolver>(comm);
  else if (solver_type == "MINRES") solver = std::make_unique<MINRESSolver>(comm);
  else if (solver_type == "GMRES")  solver = std::make_unique<GMRESSolver>(comm);
  else MFEM_ABORT("Unknown solver.");

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
// Derive from MFEM::Coefficient
class ExactTempCoeff : public Coefficient
{
public:
  double Eval(ElementTransformation &T, const IntegrationPoint &ip) override
  {
    Vector x;
    T.Transform(ip, x);
    return 270.0 + 30.0 * x[0];
  }
};