//
// Created by gangwh on 2/2/26.
//

#ifndef PCMS_MFEM_COUPLING_SCHWARTZ_COUPLING_SUPPORT_H
#define PCMS_MFEM_COUPLING_SCHWARTZ_COUPLING_SUPPORT_H

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
// -----------------------------
// Utilities
// -----------------------------
static double DefaultTolX(const Mesh &mesh)
{
  double xmin =  1e300, xmax = -1e300;
  for (int i = 0; i < mesh.GetNV(); i++)
  {
    const double *v = mesh.GetVertex(i);
    xmin = std::min(xmin, v[0]);
    xmax = std::max(xmax, v[0]);
  }
  const double Lx = xmax - xmin;
  return std::max(1e-12, 1e-10 * (std::abs(Lx) + 1.0));
}
// Extract sorted-by-y vertex samples on x=xline.
static std::vector<std::pair<double,double>>
ExtractVertexLineTrace(const ParMesh &pmesh, const ParGridFunction &T,
                       double xline, double tol)
{
  std::vector<std::pair<double,double>> trace;
  trace.reserve(pmesh.GetNV());

  for (int vi = 0; vi < pmesh.GetNV(); vi++)
  {
    const double *v = pmesh.GetVertex(vi);
    if (std::abs(v[0] - xline) <= tol)
    {
      trace.emplace_back(v[1], T(vi));
    }
  }

  std::sort(trace.begin(), trace.end(),
            [](auto &a, auto &b){ return a.first < b.first; });

  // dedup by y
  std::vector<std::pair<double,double>> uniq;
  uniq.reserve(trace.size());
  for (auto &p : trace)
  {
    if (uniq.empty() || std::abs(p.first - uniq.back().first) > 10*tol)
      uniq.push_back(p);
    else
      uniq.back().second = p.second;
  }
  return uniq;
}

// Apply a y->value trace to boundary attribute bdr_attr by setting boundary vertices (order=1).
static void ApplyBoundaryTraceByAttr(ParMesh &pmesh, ParGridFunction &gf,
                                     int bdr_attr,
                                     const std::vector<std::pair<double,double>> &trace,
                                     double tol)
{
  MFEM_VERIFY(!trace.empty(), "Empty trace passed to ApplyBoundaryTraceByAttr.");

  auto lookup = [&](double y)->double {
    auto it = std::lower_bound(trace.begin(), trace.end(), std::make_pair(y, -1e300),
                               [](auto &a, auto &b){ return a.first < b.first; });
    if (it != trace.end() && std::abs(it->first - y) <= 10*tol) return it->second;
    if (it != trace.begin())
    {
      auto it2 = std::prev(it);
      if (std::abs(it2->first - y) <= 10*tol) return it2->second;
    }
    MFEM_ABORT("Trace lookup failed: y-grid mismatch between subdomains.");
    return 0.0;
  };

  Array<int> verts;
  for (int be = 0; be < pmesh.GetNBE(); be++)
  {
    const Element *bel = pmesh.GetBdrElement(be);
    if (bel->GetAttribute() != bdr_attr) { continue; }

    pmesh.GetBdrElementVertices(be, verts);
    for (int j = 0; j < verts.Size(); j++)
    {
      const int vi = verts[j];
      const double *v = pmesh.GetVertex(vi);
      gf(vi) = lookup(v[1]);
    }
  }
}

static void ApplyBoundaryConstantByAttr(ParMesh &pmesh, ParGridFunction &gf,
                                        int bdr_attr, double value)
{
  Array<int> verts;
  for (int be = 0; be < pmesh.GetNBE(); be++)
  {
    const Element *bel = pmesh.GetBdrElement(be);
    if (bel->GetAttribute() != bdr_attr) { continue; }

    pmesh.GetBdrElementVertices(be, verts);
    for (int j = 0; j < verts.Size(); j++)
    {
      gf(verts[j]) = value;
    }
  }
}

// Diagnostics: boundary attribute min/max over boundary vertices
static void ReportBdrAttrStats(const ParMesh &pmesh, const ParGridFunction &T,
                               int bdr_attr, const char *name)
{
  Array<int> verts;
  double Tmin =  1e300, Tmax = -1e300;
  int cnt = 0;

  for (int be = 0; be < pmesh.GetNBE(); be++)
  {
    const Element *bel = pmesh.GetBdrElement(be);
    if (bel->GetAttribute() != bdr_attr) { continue; }

    pmesh.GetBdrElementVertices(be, verts);
    for (int j = 0; j < verts.Size(); j++)
    {
      const double Tv = T(verts[j]);
      Tmin = std::min(Tmin, Tv);
      Tmax = std::max(Tmax, Tv);
      cnt++;
    }
  }

  std::cout << "    " << name << " attr " << bdr_attr
            << " samples=" << cnt
            << " Tmin=" << Tmin << " Tmax=" << Tmax << "\n";
}

static void ReportTraceStats(const std::vector<std::pair<double,double>> &tr,
                             const char *name)
{
  if (tr.empty())
  {
    std::cout << "    " << name << ": EMPTY\n";
    return;
  }
  double vmin =  1e300, vmax = -1e300;
  for (auto &p : tr) { vmin = std::min(vmin, p.second); vmax = std::max(vmax, p.second); }

  std::cout << "    " << name
            << " n=" << tr.size()
            << " y=[" << tr.front().first << "," << tr.back().first << "]"
            << " val_range=[" << vmin << "," << vmax << "]"
            << " (range=" << (vmax - vmin) << ")\n";
}

static double RMSDiff(const std::vector<std::pair<double,double>> &a,
                      const std::vector<std::pair<double,double>> &b)
{
  MFEM_VERIFY(a.size() == b.size(), "Trace sizes differ.");
  double s = 0.0;
  for (size_t i = 0; i < a.size(); i++)
  {
    MFEM_VERIFY(std::abs(a[i].first - b[i].first) < 1e-10, "Trace y-grids differ.");
    const double d = a[i].second - b[i].second;
    s += d*d;
  }
  return std::sqrt(s / std::max<size_t>(1, a.size()));
}

// -----------------------------
// To help register ParaView fields
// -----------------------------
struct OutputPack
{
  ParaViewDataCollection pvd;
  ParGridFunction exact;
  ParGridFunction err;


  OutputPack(const std::string &collection, ParMesh &pm, ParFiniteElementSpace &fes)
  : pvd(collection.c_str(), &pm),
  exact(&fes),
  err(&fes)
  {
    pvd.SetDataFormat(VTKFormat::BINARY);
    pvd.SetHighOrderOutput(true);
  }
};

#endif // PCMS_MFEM_COUPLING_SCHWARTZ_COUPLING_SUPPORT_H
