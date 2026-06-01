// thermal_schwarz.cpp
//
// Steady diffusion on unit square, exact solution: T(x,y) = 270 + 30 x
//
// We solve it via overlapping Schwarz with two independent meshes:
//   ΩA = [0,0.6] x [0,1]
//   ΩB = [0.4,1] x [0,1]
//
// Iteration (Dirichlet–Dirichlet Schwarz):
//   Given gB^(k) on ΓA,right (x=0.6):
//     Solve A with T=270 on x=0 and T=gB^(k) on x=0.6
//     Extract gA^(k+1) = TA on x=0.4
//     Solve B with T=gA^(k+1) on x=0.4 and T=300 on x=1
//     Extract gB^(k+1) = TB on x=0.6
//   Stop when max(RMS(gA change), RMS(gB change)) < tol
//
// IMPORTANT:
// * This code assumes order=1 nodal H1 and that x=0.4 and x=0.6 align with vertices,
//   achieved by choosing nx multiple of 3 (default 30).
// * Boundary attribute numbering in MFEM can differ by build/version,
//   so we detect which boundary attributes correspond to x-min and x-max by geometry.
//
// Build:
//   mpicxx thermal_schwarz.cpp -I<MFEM_INC> -L<MFEM_LIB> -lmfem -lhypre -o thermal_schwarz
// Run (serial OK):
//   ./thermal_schwarz
//
// Expected:
// * boundary checks show exact constants on x-min/x-max boundaries,
// * trace values become ~constant in y,
// * errors vs exact 270+30x go to ~machine tolerance.

#include "mfem.hpp"

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
// SolveSystem (low verbosity)
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

// Detect which boundary attributes correspond to x-min and x-max boundaries.
static void FindXMinMaxBoundaryAttributes(const ParMesh &pmesh,
                                         int &attr_xmin,
                                         int &attr_xmax)
{
  double xmin =  1e300, xmax = -1e300;
  for (int vi = 0; vi < pmesh.GetNV(); vi++)
  {
    const double *v = pmesh.GetVertex(vi);
    xmin = std::min(xmin, v[0]);
    xmax = std::max(xmax, v[0]);
  }
  const double tol = std::max(1e-12, 1e-10*(std::abs(xmax-xmin)+1.0));

  std::map<int, double> sum_xavg;
  std::map<int, int> cnt;

  Array<int> verts;
  for (int be = 0; be < pmesh.GetNBE(); be++)
  {
    const Element *bel = pmesh.GetBdrElement(be);
    const int attr = bel->GetAttribute();

    pmesh.GetBdrElementVertices(be, verts);
    double xavg = 0.0;
    for (int j = 0; j < verts.Size(); j++) { xavg += pmesh.GetVertex(verts[j])[0]; }
    xavg /= verts.Size();

    sum_xavg[attr] += xavg;
    cnt[attr] += 1;
  }

  attr_xmin = -1;
  attr_xmax = -1;
  double best_min = 1e300, best_max = 1e300;

  for (auto &kv : cnt)
  {
    const int attr = kv.first;
    const double xavg = sum_xavg[attr] / kv.second;

    const double dmin = std::abs(xavg - xmin);
    const double dmax = std::abs(xavg - xmax);

    if (dmin < best_min) { best_min = dmin; attr_xmin = attr; }
    if (dmax < best_max) { best_max = dmax; attr_xmax = attr; }
  }

  MFEM_VERIFY(attr_xmin > 0 && attr_xmax > 0, "Failed to detect x-min/x-max boundary attributes.");
  MFEM_VERIFY(best_min <= 100*tol && best_max <= 100*tol,
              "Boundary attribute detection: not close enough to x-min/x-max.");
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

static std::vector<std::pair<double,double>>
RelaxTrace(const std::vector<std::pair<double,double>> &old_t,
           const std::vector<std::pair<double,double>> &new_t,
           double omega)
{
  MFEM_VERIFY(old_t.size() == new_t.size(), "Trace sizes differ.");
  std::vector<std::pair<double,double>> out = new_t;
  for (size_t i = 0; i < out.size(); i++)
  {
    out[i].second = omega*new_t[i].second + (1.0-omega)*old_t[i].second;
  }
  return out;
}

// Compare to exact T=270+30x on vertex samples
static void ErrorToExact_270_30x(const ParMesh &pmesh, const ParGridFunction &T,
                                double &rms, double &emax)
{
  double se = 0.0;
  emax = 0.0;
  const int nv = pmesh.GetNV();
  for (int vi = 0; vi < nv; vi++)
  {
    const double *v = pmesh.GetVertex(vi);
    const double Tex = 270.0 + 30.0 * v[0];
    const double e = T(vi) - Tex;
    se += e*e;
    emax = std::max(emax, std::abs(e));
  }
  rms = std::sqrt(se / std::max(1, nv));
}

// Count how many vertices lie on x=xline
static int CountLineVerts(const ParMesh &pmesh, double xline, double tol)
{
  int c = 0;
  for (int vi = 0; vi < pmesh.GetNV(); vi++)
  {
    const double *v = pmesh.GetVertex(vi);
    if (std::abs(v[0] - xline) <= tol) { c++; }
  }
  return c;
}

// Visualize the solution after every iteration
static void SaveSchwarzIteration(ParaViewDataCollection &pvdA,
                                 ParaViewDataCollection &pvdB,
                                 ParGridFunction &TA,
                                 ParGridFunction &TB,
                                 int iter,
                                 double time,
                                 bool save_exact_and_error = true,
                                 bool high_order = true)
{
  pvdA.SetCycle(iter);
  pvdA.SetTime(time);
  pvdA.SetHighOrderOutput(high_order);

  pvdB.SetCycle(iter);
  pvdB.SetTime(time);
  pvdB.SetHighOrderOutput(high_order);

  if (save_exact_and_error)
  {
    ExactTempCoeff exact;

    // Exact fields on subdomains
    ParGridFunction exactA(TA.ParFESpace());
    ParGridFunction exactB(TB.ParFESpace());
    exactA.ProjectCoefficient(exact);
    exactB.ProjectCoefficient(exact);

    // Error fields
    ParGridFunction errA(TA.ParFESpace());
    ParGridFunction errB(TB.ParFESpace());
    errA = TA; errA -= exactA;
    errB = TB; errB -= exactB;

    // Register (safe to call repeatedly; MFEM ignores duplicates by name in most builds,
    // but to be safe, you can register once outside—see notes below)
    pvdA.RegisterField("T", &TA);
    pvdA.RegisterField("T_exact", &exactA);
    pvdA.RegisterField("error", &errA);

    pvdB.RegisterField("T", &TB);
    pvdB.RegisterField("T_exact", &exactB);
    pvdB.RegisterField("error", &errB);
  }
  else
  {
    pvdA.RegisterField("T", &TA);
    pvdB.RegisterField("T", &TB);
  }

  pvdA.Save();
  pvdB.Save();
}



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

// -----------------------------
// main: Schwarz iteration + exact comparison
// -----------------------------
int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);

  {
    const int rank = Mpi::WorldRank();
    const bool is_root = (rank == 0);

    // ---- Parameters ----
    int order = 1;
    int nx = 30;         // must be multiple of 3
    int ny = 30;
    double kappa = 1.0;

    double T_left  = 270.0;
    double T_right = 300.0;

    int max_schwarz_iter = 200;
    double tol = 1e-8;
    double omega = 1.0;

    double rel_tol = 1e-12;
    int max_lin_iter = 400;

    OptionsParser args(argc, argv);
    args.AddOption(&nx, "-nx", "--nx", "Elements in x for each subdomain (must be multiple of 3).");
    args.AddOption(&ny, "-ny", "--ny", "Elements in y for each subdomain.");
    args.AddOption(&max_schwarz_iter, "-it", "--schwarz-iters", "Max Schwarz iterations.");
    args.AddOption(&tol, "-tol", "--tol", "Schwarz RMS tolerance on traces.");
    args.AddOption(&omega, "-w", "--omega", "Relaxation for trace update (0<omega<=1).");
    args.Parse();
    if (!args.Good())
    {
      if (is_root) { args.PrintUsage(std::cout); }
      MPI_Finalize();
      return 1;
    }
    if (is_root) { args.PrintOptions(std::cout); }

    MFEM_VERIFY(order == 1, "This implementation assumes order=1 nodal H1.");
    MFEM_VERIFY((nx % 3) == 0, "Choose nx multiple of 3 so x=0.4 and x=0.6 align with vertices.");

    // ---- Build subdomain meshes ----
    // ΩA: [0,0.6]x[0,1]
    Mesh smeshA = Mesh::MakeCartesian2D(nx, ny, Element::TRIANGLE, true, 0.6, 1.0);

    // ΩB: [0.4,1]x[0,1] (build [0,0.6] then shift by +0.4)
    Mesh smeshB = Mesh::MakeCartesian2D(nx, ny, Element::TRIANGLE, true, 0.6, 1.0);
    for (int i = 0; i < smeshB.GetNV(); i++) { smeshB.GetVertex(i)[0] += 0.4; }

    ParMesh pmeshA(MPI_COMM_WORLD, smeshA);
    ParMesh pmeshB(MPI_COMM_WORLD, smeshB);

    const double tolA = DefaultTolX(pmeshA);
    const double tolB = DefaultTolX(pmeshB);

    // Detect x-min/x-max boundary attributes by geometry
    int A_attr_xmin, A_attr_xmax;
    int B_attr_xmin, B_attr_xmax;
    FindXMinMaxBoundaryAttributes(pmeshA, A_attr_xmin, A_attr_xmax);
    FindXMinMaxBoundaryAttributes(pmeshB, B_attr_xmin, B_attr_xmax);

    if (is_root)
    {
      std::cout << "Detected boundary attrs:\n"
                << "  A: x-min attr=" << A_attr_xmin << "  x-max attr=" << A_attr_xmax << "  (expect x in [0,0.6])\n"
                << "  B: x-min attr=" << B_attr_xmin << "  x-max attr=" << B_attr_xmax << "  (expect x in [0.4,1])\n"
                << "  Line vertex counts:\n"
                << "    A: x=0.4 " << CountLineVerts(pmeshA, 0.4, tolA)
                << "  x=0.6 " << CountLineVerts(pmeshA, 0.6, tolA) << "\n"
                << "    B: x=0.4 " << CountLineVerts(pmeshB, 0.4, tolB)
                << "  x=0.6 " << CountLineVerts(pmeshB, 0.6, tolB) << "\n"
                << "Exact reference: T(0.4)=282, T(0.6)=288\n";
    }

    // ---- Build systems ----
    FEMSystem sysA = InitThermalSystem(&pmeshA, order, kappa);
    FEMSystem sysB = InitThermalSystem(&pmeshB, order, kappa);

    // Essential boundaries: both x-min and x-max for each subdomain
    Array<int> ess_bdrA(pmeshA.bdr_attributes.Max()); ess_bdrA = 0;
    Array<int> ess_bdrB(pmeshB.bdr_attributes.Max()); ess_bdrB = 0;

    ess_bdrA[A_attr_xmin - 1] = 1;
    ess_bdrA[A_attr_xmax - 1] = 1;
    ess_bdrB[B_attr_xmin - 1] = 1;
    ess_bdrB[B_attr_xmax - 1] = 1;

    sysA.fes->GetEssentialTrueDofs(ess_bdrA, sysA.ess_tdofs);
    sysB.fes->GetEssentialTrueDofs(ess_bdrB, sysB.ess_tdofs);

    // ---- Initialize interface trace gB on x=0.6 using A's y-grid ----
    auto ygridA = ExtractVertexLineTrace(pmeshA, *sysA.x, 0.0, tolA);
    MFEM_VERIFY(!ygridA.empty(), "Failed to get y-grid from A.");

    std::vector<std::pair<double,double>> gB_on_0p6 = ygridA;
    for (auto &p : gB_on_0p6) { p.second = 278.0; } // initial guess

    std::vector<std::pair<double,double>> gA_on_0p4_old;
    std::vector<std::pair<double,double>> gB_on_0p6_old = gB_on_0p6;

    // ---- Output ----
    OutputPack outA("schwarz_A", pmeshA, *sysA.fes);
    OutputPack outB("schwarz_B", pmeshB, *sysB.fes);

    // Register fields ONCE (pointers must remain valid)
    outA.pvd.RegisterField("T", sysA.x);
    outA.pvd.RegisterField("T_exact", &outA.exact);
    outA.pvd.RegisterField("error", &outA.err);

    outB.pvd.RegisterField("T", sysB.x);
    outB.pvd.RegisterField("T_exact", &outB.exact);
    outB.pvd.RegisterField("error", &outB.err);

    // ---- Schwarz loop ----
    for (int it = 0; it < max_schwarz_iter; it++)
    {
      // A solve: T=270 on x-min, T=gB on x-max (x=0.6)
      *sysA.x = 0.0;
      ApplyBoundaryConstantByAttr(pmeshA, *sysA.x, A_attr_xmin, T_left);
      ApplyBoundaryTraceByAttr   (pmeshA, *sysA.x, A_attr_xmax, gB_on_0p6, tolA);

      if (is_root)
      {
        std::cout << "\n[it " << it << "] Pre-solve A boundary check:\n";
        ReportBdrAttrStats(pmeshA, *sysA.x, A_attr_xmin, "A x-min");
        ReportBdrAttrStats(pmeshA, *sysA.x, A_attr_xmax, "A x-max");
        ReportTraceStats(gB_on_0p6, "gB applied on A(x=0.6)");
      }

      SolveSystem(sysA, "CG", "HypreAMG", rel_tol, max_lin_iter);

      // Extract gA_new = TA at x=0.4
      auto gA_new = ExtractVertexLineTrace(pmeshA, *sysA.x, 0.4, tolA);
      MFEM_VERIFY(!gA_new.empty(), "A has no vertices on x=0.4 (nx multiple of 3 required).");

      // B solve: T=gA on x-min (x=0.4), T=300 on x-max (x=1)
      *sysB.x = 0.0;
      ApplyBoundaryTraceByAttr   (pmeshB, *sysB.x, B_attr_xmin, gA_new, tolB);
      ApplyBoundaryConstantByAttr(pmeshB, *sysB.x, B_attr_xmax, T_right);

      if (is_root)
      {
        std::cout << "[it " << it << "] Pre-solve B boundary check:\n";
        ReportBdrAttrStats(pmeshB, *sysB.x, B_attr_xmin, "B x-min");
        ReportBdrAttrStats(pmeshB, *sysB.x, B_attr_xmax, "B x-max");
        ReportTraceStats(gA_new, "gA applied on B(x=0.4)");
      }

      SolveSystem(sysB, "CG", "HypreAMG", rel_tol, max_lin_iter);

      // Extract gB_new = TB at x=0.6
      auto gB_new = ExtractVertexLineTrace(pmeshB, *sysB.x, 0.6, tolB);
      MFEM_VERIFY(!gB_new.empty(), "B has no vertices on x=0.6 (nx multiple of 3 required).");

      // Convergence metrics BEFORE overwriting "old"
      const double rmsB = RMSDiff(gB_new, gB_on_0p6_old);
      const double rmsA = (!gA_on_0p4_old.empty()) ? RMSDiff(gA_new, gA_on_0p4_old) : rmsB;

      // Relax/update gB
      const auto gB_relaxed = RelaxTrace(gB_on_0p6_old, gB_new, omega);

      gA_on_0p4_old = gA_new;
      gB_on_0p6_old = gB_relaxed;
      gB_on_0p6     = gB_relaxed;

      // Errors to exact solution
      double rmsEA, maxEA, rmsEB, maxEB;
      ErrorToExact_270_30x(pmeshA, *sysA.x, rmsEA, maxEA);
      ErrorToExact_270_30x(pmeshB, *sysB.x, rmsEB, maxEB);

      if (is_root)
      {
        std::cout << "[it " << it << "] Trace stats after solves:\n";
        ReportTraceStats(gA_new, "gA_new (A @ x=0.4)");
        ReportTraceStats(gB_new, "gB_new (B @ x=0.6)");

        std::cout << std::setprecision(10)
                  << "[it " << it << "] rmsA(trace)=" << rmsA
                  << "  rmsB(trace)=" << rmsB
                  << "  |  A vs exact: rms=" << rmsEA << " max=" << maxEA
                  << "  B vs exact: rms=" << rmsEB << " max=" << maxEB
                  << "\n";

        // -------- SAVE FIELDS --------
        ExactTempCoeff exact;

        // update exact
        outA.exact.ProjectCoefficient(exact);
        outB.exact.ProjectCoefficient(exact);

        // update error = T - T_exact
        outA.err = *sysA.x; outA.err -= outA.exact;
        outB.err = *sysB.x; outB.err -= outB.exact;

        // time-series metadata
        outA.pvd.SetCycle(it);
        outA.pvd.SetTime((double)it);
        outB.pvd.SetCycle(it);
        outB.pvd.SetTime((double)it);

        // write
        outA.pvd.Save();
        outB.pvd.Save();
      }

      if (std::max(rmsA, rmsB) < tol)
      {
        if (is_root) { std::cout << "Converged (max RMS trace < tol).\n"; }
        break;
      }
    }

    if (is_root)
    {
      std::cout << "\nFinal relaxed gB(x=0.6):\n";
      ReportTraceStats(gB_on_0p6_old, "gB(x=0.6)");
    }

    // ---- Clean up ----
    delete sysA.x; delete sysA.a; delete sysA.b; delete sysA.fes; delete sysA.fec;
    delete sysB.x; delete sysB.a; delete sysB.b; delete sysB.fes; delete sysB.fec;
  }

  MPI_Finalize();
  return 0;
}