//
// Created by gangwh on 4/22/26.
//
#include "Schwarz_Sundial_Coupling.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <stdexcept>
#include <utility>

namespace schwarz {

using namespace mfem;

namespace {

SUNErrCode SchwarzStepper_Evolve(SUNStepper stepper,
                                 sunrealtype tout,
                                 N_Vector vret,
                                 sunrealtype* tret)
{
  void* content_void = nullptr;
  if (SUNStepper_GetContent(stepper, &content_void) != SUN_SUCCESS)
  {
    return SUN_ERR_EXT_FAIL;
  }

  auto* C = static_cast<SchwarzStepperContent*>(content_void);

  try
  {
    Trace gA = C->gA_meta;
    Trace gB = C->gB_meta;
    UnpackState(vret, gA, gB);

    while (C->tcur < tout && C->tcur < C->cfg.pseudo_tstop)
    {
      Trace gA_new = C->gA_meta;
      Trace gB_new = C->gB_meta;

      const int flag = SchwarzSweep(C, gA, gB, gA_new, gB_new);
      if (flag != 0)
      {
        SUNStepper_SetLastFlag(stepper, flag);
        return SUN_ERR_EXT_FAIL;
      }

      gA = std::move(gA_new);
      gB = std::move(gB_new);
      C->tcur += C->cfg.pseudo_dt;
      C->nsteps += 1;
    }

    PackState(gA, gB, vret);
    if (tret) { *tret = C->tcur; }

    SUNStepper_SetLastFlag(stepper, 0);
    return SUN_SUCCESS;
  }
  catch (...)
  {
    SUNStepper_SetLastFlag(stepper, -1);
    return SUN_ERR_EXT_FAIL;
  }
}

SUNErrCode SchwarzStepper_Reset(SUNStepper stepper,
                                sunrealtype tR,
                                N_Vector /* vR */)
{
  void* content_void = nullptr;
  if (SUNStepper_GetContent(stepper, &content_void) != SUN_SUCCESS)
  {
    return SUN_ERR_EXT_FAIL;
  }

  auto* C = static_cast<SchwarzStepperContent*>(content_void);
  C->tcur = tR;
  C->nsteps = 0;
  return SUN_SUCCESS;
}

SUNErrCode SchwarzStepper_ReInit(SUNStepper stepper,
                                 sunrealtype tR,
                                 N_Vector vR)
{
  return SchwarzStepper_Reset(stepper, tR, vR);
}

SUNErrCode SchwarzStepper_SetStopTime(SUNStepper stepper, sunrealtype tstop)
{
  void* content_void = nullptr;
  if (SUNStepper_GetContent(stepper, &content_void) != SUN_SUCCESS)
  {
    return SUN_ERR_EXT_FAIL;
  }

  auto* C = static_cast<SchwarzStepperContent*>(content_void);
  C->cfg.pseudo_tstop = tstop;
  return SUN_SUCCESS;
}

SUNErrCode SchwarzStepper_GetNumSteps(SUNStepper stepper, suncountertype* nsteps)
{
  void* content_void = nullptr;
  if (SUNStepper_GetContent(stepper, &content_void) != SUN_SUCCESS)
  {
    return SUN_ERR_EXT_FAIL;
  }

  auto* C = static_cast<SchwarzStepperContent*>(content_void);
  *nsteps = C->nsteps;
  return SUN_SUCCESS;
}

static SUNErrCode SchwarzStepper_Destroy(SUNStepper stepper)
{
  void* content_void = nullptr;
  if (SUNStepper_GetContent(stepper, &content_void) != SUN_SUCCESS)
    return SUN_ERR_EXT_FAIL;

  auto* C = static_cast<SchwarzStepperContent*>(content_void);
  if (!C) return SUN_SUCCESS;

  delete C->sysA.x;
  delete C->sysA.a;
  delete C->sysA.b;
  delete C->sysA.fes;
  delete C->sysA.fec;

  delete C->sysB.x;
  delete C->sysB.a;
  delete C->sysB.b;
  delete C->sysB.fes;
  delete C->sysB.fec;

  delete C;

  // optional but clean
  SUNStepper_SetContent(stepper, nullptr);

  return SUN_SUCCESS;
}
} // namespace

FEMSystem InitThermalSystem(ParMesh* pmesh, int order, double kappa_val)
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

void DestroyFEMSystem(FEMSystem& sys)
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

double DefaultTolX(const Mesh& mesh)
{
  double xmin = 1e300;
  double xmax = -1e300;
  for (int i = 0; i < mesh.GetNV(); i++)
  {
    const double* v = mesh.GetVertex(i);
    xmin = std::min(xmin, v[0]);
    xmax = std::max(xmax, v[0]);
  }

  const double Lx = xmax - xmin;
  return std::max(1e-12, 1e-10 * (std::abs(Lx) + 1.0));
}

void FindXMinMaxBoundaryAttributes(const ParMesh& pmesh,
                                   int& attr_xmin,
                                   int& attr_xmax)
{
  double xmin = 1e300;
  double xmax = -1e300;
  for (int vi = 0; vi < pmesh.GetNV(); vi++)
  {
    const double* v = pmesh.GetVertex(vi);
    xmin = std::min(xmin, v[0]);
    xmax = std::max(xmax, v[0]);
  }

  const double tol = std::max(1e-12, 1e-10 * (std::abs(xmax - xmin) + 1.0));

  std::map<int, double> sum_xavg;
  std::map<int, int> cnt;

  Array<int> verts;
  for (int be = 0; be < pmesh.GetNBE(); be++)
  {
    const Element* bel = pmesh.GetBdrElement(be);
    const int attr = bel->GetAttribute();

    pmesh.GetBdrElementVertices(be, verts);
    double xavg = 0.0;
    for (int j = 0; j < verts.Size(); j++)
    {
      xavg += pmesh.GetVertex(verts[j])[0];
    }
    xavg /= verts.Size();

    sum_xavg[attr] += xavg;
    cnt[attr] += 1;
  }

  attr_xmin = -1;
  attr_xmax = -1;
  double best_min = 1e300;
  double best_max = 1e300;

  for (const auto& kv : cnt)
  {
    const int attr = kv.first;
    const double xavg = sum_xavg[attr] / kv.second;

    const double dmin = std::abs(xavg - xmin);
    const double dmax = std::abs(xavg - xmax);

    if (dmin < best_min)
    {
      best_min = dmin;
      attr_xmin = attr;
    }
    if (dmax < best_max)
    {
      best_max = dmax;
      attr_xmax = attr;
    }
  }

  MFEM_VERIFY(attr_xmin > 0 && attr_xmax > 0,
              "Failed to detect x-min/x-max boundary attributes.");
  MFEM_VERIFY(best_min <= 100 * tol && best_max <= 100 * tol,
              "Boundary attribute detection not close enough to x-min/x-max.");
}

Trace PairTraceToTrace(const std::vector<std::pair<double, double>>& in)
{
  Trace out;
  out.y.reserve(in.size());
  out.val.reserve(in.size());

  for (const auto& p : in)
  {
    out.y.push_back(p.first);
    out.val.push_back(p.second);
  }

  return out;
}

std::vector<std::pair<double, double>> TraceToPairTrace(const Trace& t)
{
  std::vector<std::pair<double, double>> out(t.val.size());
  for (int i = 0; i < static_cast<int>(t.val.size()); i++)
  {
    out[i] = {t.y[i], t.val[i]};
  }
  return out;
}

Trace BlendTrace(const Trace& gA, const Trace& gB, double omega)
{
  if (gA.Size() != gB.Size())
  {
    throw std::runtime_error("Trace size mismatch in BlendTrace");
  }

  Trace G = gA;
  for (int i = 0; i < G.Size(); i++)
  {
    if (std::abs(gA.y[i] - gB.y[i]) > 1e-12)
    {
      throw std::runtime_error("Trace y-grid mismatch in BlendTrace");
    }
    G.val[i] = omega * gA.val[i] + (1.0 - omega) * gB.val[i];
  }

  return G;
}

std::vector<std::pair<double, double>>
ExtractVertexLineTrace(const ParMesh& pmesh,
                       const ParGridFunction& T,
                       double xline,
                       double tol)
{
  std::vector<std::pair<double, double>> trace;
  trace.reserve(pmesh.GetNV());

  for (int vi = 0; vi < pmesh.GetNV(); vi++)
  {
    const double* v = pmesh.GetVertex(vi);
    if (std::abs(v[0] - xline) <= tol)
    {
      trace.emplace_back(v[1], T(vi));
    }
  }

  std::sort(trace.begin(), trace.end(),
            [](const auto& a, const auto& b) { return a.first < b.first; });

  std::vector<std::pair<double, double>> uniq;
  uniq.reserve(trace.size());
  for (const auto& p : trace)
  {
    if (uniq.empty() || std::abs(p.first - uniq.back().first) > 10 * tol)
    {
      uniq.push_back(p);
    }
    else
    {
      uniq.back().second = p.second;
    }
  }

  return uniq;
}

void ApplyBoundaryTraceByAttr(ParMesh& pmesh,
                              ParGridFunction& gf,
                              int bdr_attr,
                              const std::vector<std::pair<double, double>>& trace,
                              double tol)
{
  MFEM_VERIFY(!trace.empty(), "Empty trace passed to ApplyBoundaryTraceByAttr.");

  auto lookup = [&](double y) -> double {
    auto it = std::lower_bound(trace.begin(), trace.end(), std::make_pair(y, -1e300),
                               [](const auto& a, const auto& b) {
                                 return a.first < b.first;
                               });
    if (it != trace.end() && std::abs(it->first - y) <= 10 * tol)
    {
      return it->second;
    }
    if (it != trace.begin())
    {
      auto it2 = std::prev(it);
      if (std::abs(it2->first - y) <= 10 * tol)
      {
        return it2->second;
      }
    }
    MFEM_ABORT("Trace lookup failed: y-grid mismatch between subdomains.");
    return 0.0;
  };

  Array<int> verts;
  for (int be = 0; be < pmesh.GetNBE(); be++)
  {
    const Element* bel = pmesh.GetBdrElement(be);
    if (bel->GetAttribute() != bdr_attr) { continue; }

    pmesh.GetBdrElementVertices(be, verts);
    for (int j = 0; j < verts.Size(); j++)
    {
      const int vi = verts[j];
      const double* v = pmesh.GetVertex(vi);
      gf(vi) = lookup(v[1]);
    }
  }
}

void ApplyBoundaryConstantByAttr(ParMesh& pmesh,
                                 ParGridFunction& gf,
                                 int bdr_attr,
                                 double value)
{
  Array<int> verts;
  for (int be = 0; be < pmesh.GetNBE(); be++)
  {
    const Element* bel = pmesh.GetBdrElement(be);
    if (bel->GetAttribute() != bdr_attr) { continue; }

    pmesh.GetBdrElementVertices(be, verts);
    for (int j = 0; j < verts.Size(); j++)
    {
      gf(verts[j]) = value;
    }
  }
}

void PackState(const Trace& gA, const Trace& gB, N_Vector nv)
{
  sunrealtype* data = N_VGetArrayPointer(nv);
  const sunindextype n = N_VGetLength(nv);
  const sunindextype m = n / 2;

  if (n % 2 != 0)
  {
    throw std::runtime_error("State vector length must be even in PackState");
  }
  if (gA.Size() != m || gB.Size() != m)
  {
    throw std::runtime_error("Trace size/state length mismatch in PackState");
  }

  for (sunindextype i = 0; i < m; i++)
  {
    data[i] = gA.val[i];
    data[m + i] = gB.val[i];
  }
}

void UnpackState(N_Vector nv, Trace& gA, Trace& gB)
{
  sunrealtype* data = N_VGetArrayPointer(nv);
  const sunindextype n = N_VGetLength(nv);
  if (n % 2 != 0)
  {
    throw std::runtime_error("State vector length must be even in UnpackState");
  }

  const sunindextype m = n / 2;
  if (gA.Size() != m || gB.Size() != m)
  {
    throw std::runtime_error("Trace size/state length mismatch in UnpackState");
  }

  for (sunindextype i = 0; i < m; i++)
  {
    gA.val[i] = data[i];
    gB.val[i] = data[m + i];
  }
}

double SolveSystem(FEMSystem& sys,
                   const std::string& solver_type,
                   const std::string& prec_type,
                   double rel_tol,
                   int max_iter)
{
  OperatorPtr A;
  HypreParVector X, B;

  sys.a->FormLinearSystem(sys.ess_tdofs, *sys.x, *sys.b, A, X, B);

  auto* A_hypre = A.As<HypreParMatrix>();
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

  if (solver_type == "CG")
  {
    solver = std::make_unique<CGSolver>(comm);
  }
  else if (solver_type == "MINRES")
  {
    solver = std::make_unique<MINRESSolver>(comm);
  }
  else if (solver_type == "GMRES")
  {
    solver = std::make_unique<GMRESSolver>(comm);
  }
  else
  {
    MFEM_ABORT("Unknown solver.");
  }

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

int SchwarzSweep(SchwarzStepperContent* C,
                 const Trace& gA_old,
                 const Trace& gB_old,
                 Trace& gA_new,
                 Trace& gB_new)
{
  Trace G = BlendTrace(gA_old, gB_old, C->cfg.omega);

  *C->sysA.x = 0.0;
  ApplyBoundaryConstantByAttr(*C->sysA.pmesh, *C->sysA.x, C->A_attr_xmin,
                              C->cfg.T_left);
  ApplyBoundaryTraceByAttr(*C->sysA.pmesh, *C->sysA.x, C->A_attr_xmax,
                           TraceToPairTrace(G), C->tolA);

  SolveSystem(C->sysA, C->cfg.solver_type, C->cfg.prec_type,
              C->cfg.rel_tol, C->cfg.max_lin_iter);

  gA_new = PairTraceToTrace(
      ExtractVertexLineTrace(*C->sysA.pmesh, *C->sysA.x,
                             C->cfg.x_A_extract, C->tolA));
  if (gA_new.Empty())
  {
    throw std::runtime_error("Extracted empty gA trace in SchwarzSweep");
  }

  *C->sysB.x = 0.0;
  ApplyBoundaryTraceByAttr(*C->sysB.pmesh, *C->sysB.x, C->B_attr_xmin,
                           TraceToPairTrace(gA_new), C->tolB);
  ApplyBoundaryConstantByAttr(*C->sysB.pmesh, *C->sysB.x, C->B_attr_xmax,
                              C->cfg.T_right);

  SolveSystem(C->sysB, C->cfg.solver_type, C->cfg.prec_type,
              C->cfg.rel_tol, C->cfg.max_lin_iter);

  gB_new = PairTraceToTrace(
      ExtractVertexLineTrace(*C->sysB.pmesh, *C->sysB.x,
                             C->cfg.x_B_extract, C->tolB));
  if (gB_new.Empty())
  {
    throw std::runtime_error("Extracted empty gB trace in SchwarzSweep");
  }

  if (gA_new.Size() != C->gA_meta.Size() || gB_new.Size() != C->gB_meta.Size())
  {
    throw std::runtime_error("Trace size changed across Schwarz sweeps");
  }

  return 0;
}

std::unique_ptr<SchwarzStepperContent>
BuildDefaultSchwarzContent(MPI_Comm comm, const SchwarzConfig& cfg)
{
  MFEM_VERIFY(cfg.order == 1,
              "This template assumes order=1 nodal H1 elements.");
  MFEM_VERIFY((cfg.nx % 3) == 0,
              "Choose nx multiple of 3 so x=0.4 and x=0.6 align with vertices.");

  auto content = std::make_unique<SchwarzStepperContent>();
  content->cfg = cfg;
  content->tcur = cfg.pseudo_t0;

  Mesh smeshA = Mesh::MakeCartesian2D(cfg.nx, cfg.ny, Element::TRIANGLE,
                                      true, 0.6, 1.0);
  Mesh smeshB = Mesh::MakeCartesian2D(cfg.nx, cfg.ny, Element::TRIANGLE,
                                      true, 0.6, 1.0);
  for (int i = 0; i < smeshB.GetNV(); i++)
  {
    smeshB.GetVertex(i)[0] += 0.4;
  }

  auto* pmeshA = new ParMesh(comm, smeshA);
  auto* pmeshB = new ParMesh(comm, smeshB);

  content->sysA = InitThermalSystem(pmeshA, cfg.order, cfg.kappa);
  content->sysB = InitThermalSystem(pmeshB, cfg.order, cfg.kappa);

  content->tolA = DefaultTolX(*pmeshA);
  content->tolB = DefaultTolX(*pmeshB);

  FindXMinMaxBoundaryAttributes(*pmeshA, content->A_attr_xmin, content->A_attr_xmax);
  FindXMinMaxBoundaryAttributes(*pmeshB, content->B_attr_xmin, content->B_attr_xmax);

  Array<int> ess_bdrA(pmeshA->bdr_attributes.Max());
  ess_bdrA = 0;
  ess_bdrA[content->A_attr_xmin - 1] = 1;
  ess_bdrA[content->A_attr_xmax - 1] = 1;
  content->sysA.fes->GetEssentialTrueDofs(ess_bdrA, content->sysA.ess_tdofs);

  Array<int> ess_bdrB(pmeshB->bdr_attributes.Max());
  ess_bdrB = 0;
  ess_bdrB[content->B_attr_xmin - 1] = 1;
  ess_bdrB[content->B_attr_xmax - 1] = 1;
  content->sysB.fes->GetEssentialTrueDofs(ess_bdrB, content->sysB.ess_tdofs);

  content->gA_meta = PairTraceToTrace(
      ExtractVertexLineTrace(*pmeshA, *content->sysA.x, cfg.x_A_extract, content->tolA));
  content->gB_meta = PairTraceToTrace(
      ExtractVertexLineTrace(*pmeshB, *content->sysB.x, cfg.x_B_extract, content->tolB));

  MFEM_VERIFY(!content->gA_meta.Empty(), "Failed to build A trace metadata.");
  MFEM_VERIFY(!content->gB_meta.Empty(), "Failed to build B trace metadata.");
  MFEM_VERIFY(content->gA_meta.Size() == content->gB_meta.Size(),
              "A and B trace sizes differ; expected a consistent interface grid.");

  for (int i = 0; i < content->gA_meta.Size(); i++)
  {
    MFEM_VERIFY(std::abs(content->gA_meta.y[i] - content->gB_meta.y[i]) < 1e-10,
                "A and B trace y-grids differ.");
  }

  return content;
}

SUNErrCode CreateSchwarzSUNStepper(SUNContext sunctx,
                                   SchwarzStepperContent* content,
                                   SUNStepper* stepper)
{
  SUNErrCode err = SUNStepper_Create(sunctx, stepper);
  if (err != SUN_SUCCESS) { return err; }

  err = SUNStepper_SetContent(*stepper, content);
  if (err != SUN_SUCCESS) { return err; }

  err = SUNStepper_SetEvolveFn(*stepper, SchwarzStepper_Evolve);
  if (err != SUN_SUCCESS) { return err; }

  err = SUNStepper_SetResetFn(*stepper, SchwarzStepper_Reset);
  if (err != SUN_SUCCESS) { return err; }

  err = SUNStepper_SetReInitFn(*stepper, SchwarzStepper_ReInit);
  if (err != SUN_SUCCESS) { return err; }

  err = SUNStepper_SetStopTimeFn(*stepper, SchwarzStepper_SetStopTime);
  if (err != SUN_SUCCESS) { return err; }

  err = SUNStepper_SetGetNumStepsFn(*stepper, SchwarzStepper_GetNumSteps);
  if (err != SUN_SUCCESS) { return err; }

  err = SUNStepper_SetDestroyFn(*stepper, SchwarzStepper_Destroy);
  if (err != SUN_SUCCESS) { return err; }

  return SUN_SUCCESS;
}

} // namespace schwarz