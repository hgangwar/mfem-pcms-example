#include "include/coupling_support.h"

namespace support
{

FEMSystem Init_FEMSystem(mfem::ParMesh* pmesh, int order, double kappa_val)
{
  FEMSystem sys;
  sys.pmesh = pmesh;

  const int dim = pmesh->Dimension();
  sys.fec = new mfem::H1_FECollection(order, dim);
  sys.fes = new mfem::ParFiniteElementSpace(pmesh, sys.fec);

  mfem::ConstantCoefficient zero(0.0);
  sys.b = new mfem::ParLinearForm(sys.fes);
  sys.b->AddDomainIntegrator(new mfem::DomainLFIntegrator(zero));
  sys.b->Assemble();

  mfem::ConstantCoefficient kappa(kappa_val);
  sys.a = new mfem::ParBilinearForm(sys.fes);
  sys.a->AddDomainIntegrator(new mfem::DiffusionIntegrator(kappa));
  sys.a->Assemble();
  sys.a->Finalize();

  sys.x = new mfem::ParGridFunction(sys.fes);
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

long SolveSystem(FEMSystem& sys, const std::string& solver_type,
                 const std::string& prec_type, double rel_tol, int max_iter)
{
  mfem::OperatorPtr A;
  mfem::HypreParVector X, B;

  sys.a->FormLinearSystem(sys.ess_tdofs, *sys.x, *sys.b, A, X, B);

  auto* A_hypre = A.As<mfem::HypreParMatrix>();
  MFEM_VERIFY(A_hypre, "FormLinearSystem did not produce HypreParMatrix.");

  std::unique_ptr<mfem::Solver> prec;
  if (prec_type == "HypreAMG") {
    auto amg = std::make_unique<mfem::HypreBoomerAMG>(*A_hypre);
    amg->SetPrintLevel(0);
    prec = std::move(amg);
  } else if (prec_type == "Jacobi") {
    auto sm = std::make_unique<mfem::HypreSmoother>(*A_hypre);
    sm->SetType(mfem::HypreSmoother::Jacobi);
    prec = std::move(sm);
  } else {
    MFEM_ABORT("Unknown preconditioner.");
  }

  std::unique_ptr<mfem::IterativeSolver> solver;
  MPI_Comm comm = sys.fes->GetParMesh()->GetComm();

  if (solver_type == "CG")
    solver = std::make_unique<mfem::CGSolver>(comm);
  else if (solver_type == "MINRES")
    solver = std::make_unique<mfem::MINRESSolver>(comm);
  else if (solver_type == "GMRES")
    solver = std::make_unique<mfem::GMRESSolver>(comm);
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

double ExactTempCoeff::Eval(mfem::ElementTransformation& T,
                            const mfem::IntegrationPoint& ip)
{
  mfem::Vector x;
  T.Transform(ip, x);
  return 270.0 + 30.0 * x[0];
}

double DefaultTolX(const mfem::Mesh& mesh)
{
  double xmin = 1e300, xmax = -1e300;
  for (int i = 0; i < mesh.GetNV(); i++) {
    const double* v = mesh.GetVertex(i);
    xmin = std::min(xmin, v[0]);
    xmax = std::max(xmax, v[0]);
  }
  const double Lx = xmax - xmin;
  return std::max(1e-12, 1e-10 * (std::abs(Lx) + 1.0));
}

void shift_meshX(Omega_h::Mesh& mesh, double dx)
{
  auto coords_A = Omega_h::Read<Omega_h::Real>(mesh.coords());
  auto coords_B = Omega_h::Write<Omega_h::Real>(mesh.coords().size());
  const auto nverts = mesh.nverts();
  const auto dim = mesh.dim();
  OMEGA_H_CHECK(coords_A.size() == dim * nverts);
  Omega_h::parallel_for(
    nverts, OMEGA_H_LAMBDA(const Omega_h::LO vert_id) {
      coords_B[vert_id * dim + 0] = coords_A[vert_id * dim] + dx;
    });
  mesh.set_coords(coords_B);
}

Omega_h::Write<Omega_h::I8> create_mask(Omega_h::Mesh& mesh,
                                        const char* tag_name, int tag_value)
{
  Omega_h::Write<Omega_h::I8> mask(mesh.nents(0), 0);

  const int dim = mesh.dim();
  auto tag = mesh.get_array<int>(dim, tag_name);
  std::cout << "Element tag size: " << tag.size() << std::endl;

  const auto elem2verts = mesh.ask_elem_verts();
  int unique_vertices_filtered = 0;

  if (dim == 2) {
    for (Omega_h::LO e = 0; e < mesh.nelems(); ++e) {
      if (tag[e] != tag_value) {
        continue;
      }

      const auto tri_verts = Omega_h::gather_verts<3>(elem2verts, e);
      for (int j = 0; j < 3; ++j) {
        const auto v = tri_verts[j];
        if (mask[v] == 0) {
          mask[v] = 1;
          ++unique_vertices_filtered;
        }
      }
    }
  } else if (dim == 3) {
    for (Omega_h::LO e = 0; e < mesh.nelems(); ++e) {
      if (tag[e] != tag_value) {
        continue;
      }

      const auto tet_verts = Omega_h::gather_verts<4>(elem2verts, e);
      for (int j = 0; j < 4; ++j) {
        const auto v = tet_verts[j];
        if (mask[v] == 0) {
          mask[v] = 1;
          ++unique_vertices_filtered;
        }
      }
    }
  } else {
    std::cerr << "Unsupported mesh dimension: " << dim << std::endl;
  }

  std::cout << "Number of unique vertices filtered: "
            << unique_vertices_filtered << std::endl;

  return mask;
}

double RMSDiff(const std::vector<std::pair<double, double>>& a,
               const std::vector<std::pair<double, double>>& b)
{
  MFEM_VERIFY(a.size() == b.size(), "Trace sizes differ.");
  double s = 0.0;
  for (size_t i = 0; i < a.size(); i++) {
    MFEM_VERIFY(std::abs(a[i].first - b[i].first) < 1e-10,
                "Trace y-grids differ.");
    const double d = a[i].second - b[i].second;
    s += d * d;
  }
  return std::sqrt(s / std::max<size_t>(1, a.size()));
}

long ComputeRMS(const Omega_h::Read<dtype>& a, const Omega_h::Read<dtype>& b)
{
  const int n = a.size();

  if (b.size() != n) {
    throw std::runtime_error("ComputeRMS: size mismatch");
  }

  if (n == 0) {
    return 0.0;
  }

  Omega_h::HostRead<dtype> ha(a);
  Omega_h::HostRead<dtype> hb(b);

  double sum_sq = 0.0;

  for (int i = 0; i < n; ++i) {
    double diff = ha[i] - hb[i];
    sum_sq += diff * diff;
  }

  return std::sqrt(sum_sq / static_cast<dtype>(n));
}

std::vector<std::pair<double, double>> ExtractVertexLineTrace(
  const Omega_h::Mesh& mesh, Omega_h::Read<Omega_h::Real> field_v, double xline,
  double tol)
{
  std::vector<std::pair<double, double>> trace;
  trace.reserve(static_cast<std::size_t>(mesh.nverts()));

  auto const coords = mesh.coords();
  int const dim = mesh.dim();
  OMEGA_H_CHECK(dim >= 2);

  for (int vi = 0; vi < mesh.nverts(); ++vi) {
    double const x = coords[vi * dim + 0];
    double const y = coords[vi * dim + 1];

    if (std::abs(x - xline) <= tol) {
      trace.emplace_back(y, static_cast<double>(field_v[vi]));
    }
  }

  std::sort(trace.begin(), trace.end(),
            [](auto const& a, auto const& b) { return a.first < b.first; });

  std::vector<std::pair<double, double>> uniq;
  uniq.reserve(trace.size());
  for (auto const& p : trace) {
    if (uniq.empty() || std::abs(p.first - uniq.back().first) > 10 * tol) {
      uniq.push_back(p);
    } else {
      uniq.back().second = p.second;
    }
  }
  return uniq;
}

void FillTagOnXLineFromTrace(
  Omega_h::Mesh& mesh, const std::vector<std::pair<double, double>>& trace,
  double x_line, double tol, const char* tag_name)
{
  if (trace.empty()) {
    return;
  }
  if (tol <= 0.0) {
    throw std::runtime_error("y_tol must be > 0");
  }

  if (!mesh.has_tag(0, tag_name)) {
    throw std::runtime_error(std::string("Mesh missing vertex tag: ") +
                             tag_name);
  }

  const int dim = mesh.dim();
  if (dim < 1) {
    throw std::runtime_error("Mesh dim invalid");
  }

  const int nverts = mesh.nverts();
  if (nverts == 0) {
    return;
  }

  std::unordered_map<long long, double> ybin_to_val;
  ybin_to_val.reserve(trace.size() * 2);

  auto ykey = [&](double y) -> long long { return llround(y / tol); };

  for (const auto& p : trace) {
    ybin_to_val[ykey(p.first)] = p.second;
  }

  Omega_h::HostRead<Omega_h::Real> hcoords(mesh.coords());

  Omega_h::Reals temp_in = mesh.get_array<Omega_h::Real>(0, tag_name);
  Omega_h::HostRead<Omega_h::Real> htemp_in(temp_in);

  Omega_h::Write<Omega_h::Real> temp_out_w(nverts);
  Omega_h::HostWrite<Omega_h::Real> htemp_out(temp_out_w);

  for (int v = 0; v < nverts; ++v) {
    htemp_out[v] = htemp_in[v];
  }

  for (int v = 0; v < nverts; ++v) {
    const double x = hcoords[v * dim + 0];
    if (std::abs(x - x_line) > tol) {
      continue;
    }

    const double y = (dim > 1) ? hcoords[v * dim + 1] : 0.0;
    const auto it = ybin_to_val.find(ykey(y));
    if (it != ybin_to_val.end()) {
      htemp_out[v] = it->second;
    }
  }

  mesh.set_tag(0, tag_name, Omega_h::Reals(temp_out_w));
}

void ApplyBoundaryTraceByAttr(
  mfem::ParMesh& pmesh, mfem::ParGridFunction& gf, int bdr_attr,
  const std::vector<std::pair<double, double>>& trace, double tol)
{
  MFEM_VERIFY(!trace.empty(),
              "Empty trace passed to ApplyBoundaryTraceByAttr.");

  auto lookup = [&](double y) -> double {
    auto it =
      std::lower_bound(trace.begin(), trace.end(), std::make_pair(y, -1e300),
                       [](auto& a, auto& b) { return a.first < b.first; });
    if (it != trace.end() && std::abs(it->first - y) <= 10 * tol) {
      return it->second;
    }
    if (it != trace.begin()) {
      auto it2 = std::prev(it);
      if (std::abs(it2->first - y) <= 10 * tol) {
        return it2->second;
      }
    }
    MFEM_ABORT("Trace lookup failed: y-grid mismatch between subdomains.");
    return 0.0;
  };

  mfem::Array<int> verts;
  for (int be = 0; be < pmesh.GetNBE(); be++) {
    const mfem::Element* bel = pmesh.GetBdrElement(be);
    if (bel->GetAttribute() != bdr_attr) {
      continue;
    }

    pmesh.GetBdrElementVertices(be, verts);
    for (int j = 0; j < verts.Size(); j++) {
      const int vi = verts[j];
      const double* v = pmesh.GetVertex(vi);
      gf(vi) = lookup(v[1]);
    }
  }
}

void ApplyBoundaryConstantByAttr(mfem::ParMesh& pmesh,
                                 mfem::ParGridFunction& gf, int bdr_attr,
                                 double value)
{
  mfem::Array<int> verts;
  for (int be = 0; be < pmesh.GetNBE(); be++) {
    const mfem::Element* bel = pmesh.GetBdrElement(be);
    if (bel->GetAttribute() != bdr_attr) {
      continue;
    }

    pmesh.GetBdrElementVertices(be, verts);
    for (int j = 0; j < verts.Size(); j++) {
      gf(verts[j]) = value;
    }
  }
}

void ReportBdrAttrStats(const mfem::ParMesh& pmesh,
                        const mfem::ParGridFunction& T, int bdr_attr,
                        const char* name)
{
  mfem::Array<int> verts;
  double Tmin = 1e300, Tmax = -1e300;
  int cnt = 0;

  for (int be = 0; be < pmesh.GetNBE(); be++) {
    const mfem::Element* bel = pmesh.GetBdrElement(be);
    if (bel->GetAttribute() != bdr_attr) {
      continue;
    }

    pmesh.GetBdrElementVertices(be, verts);
    for (int j = 0; j < verts.Size(); j++) {
      const double Tv = T(verts[j]);
      Tmin = std::min(Tmin, Tv);
      Tmax = std::max(Tmax, Tv);
      cnt++;
    }
  }

  std::cout << "    " << name << " attr " << bdr_attr << " samples=" << cnt
            << " Tmin=" << Tmin << " Tmax=" << Tmax << "\n";
}

void ReportTraceStats(const std::vector<std::pair<double, double>>& tr,
                      const char* name)
{
  if (tr.empty()) {
    std::cout << "    " << name << ": EMPTY\n";
    return;
  }
  double vmin = 1e300, vmax = -1e300;
  for (auto& p : tr) {
    vmin = std::min(vmin, p.second);
    vmax = std::max(vmax, p.second);
  }

  std::cout << "    " << name << " n=" << tr.size() << " y=["
            << tr.front().first << "," << tr.back().first << "]"
            << " val_range=[" << vmin << "," << vmax << "]"
            << " (range=" << (vmax - vmin) << ")\n";
}

OutputPack::OutputPack(const std::string& collection, mfem::ParMesh& pm,
                       mfem::ParFiniteElementSpace& fes)
  : pvd(collection.c_str(), &pm), exact(&fes), err(&fes)
{
  pvd.SetDataFormat(mfem::VTKFormat::BINARY);
  pvd.SetHighOrderOutput(true);
}

std::vector<std::pair<double, double>> RelaxTrace(
  const std::vector<std::pair<double, double>>& old_t,
  const std::vector<std::pair<double, double>>& new_t, double omega)
{
  MFEM_VERIFY(old_t.size() == new_t.size(), "Trace sizes differ.");
  std::vector<std::pair<double, double>> out = new_t;
  for (size_t i = 0; i < out.size(); i++) {
    out[i].second = omega * new_t[i].second + (1.0 - omega) * old_t[i].second;
  }
  return out;
}

void SaveFields(OutputPack& out, const FEMSystem& sys, int it)
{
  ExactTempCoeff exact;

  out.exact.ProjectCoefficient(exact);

  out.err = *sys.x;
  out.err -= out.exact;

  out.pvd.SetCycle(it);
  out.pvd.SetTime(static_cast<double>(it));

  out.pvd.Save();
}

void write_oh_mesh(Omega_h::Mesh& mesh, const std::string& path)
{
  Omega_h::binary::write(path, &mesh);
}

mfem::Mesh read_mfem_mesh(const std::string& path)
{
  return mfem::Mesh(path.c_str(), 1, 1);
}

void reset_mfem_attributes(mfem::Mesh& mesh, int attr)
{
  for (int e = 0; e < mesh.GetNE(); ++e) {
    mesh.SetAttribute(e, attr);
  }

  mesh.SetAttributes();
}

void remove_oh_tag(Omega_h::Mesh& mesh, const std::string& name)
{
  int dim = mesh.dim();

  if (mesh.has_tag(dim, name)) {
    mesh.remove_tag(dim, name);
  }
}

void SaveParaview(mfem::ParMesh& pmesh, mfem::ParGridFunction& x,
                  const std::string& collection, const std::string& field_name,
                  int cycle, double time)
{
  mfem::ParaViewDataCollection pvdc(collection.c_str(), &pmesh);
  pvdc.SetPrefixPath("paraview");
  pvdc.SetCycle(cycle);
  pvdc.SetTime(time);
  pvdc.SetHighOrderOutput(true);
  pvdc.SetDataFormat(mfem::VTKFormat::BINARY);
  pvdc.RegisterField(field_name.c_str(), &x);
  pvdc.Save();
}
double ComputeAbsoluteError(const Omega_h::Mesh& mesh)
{
  auto temp = mesh.get_array<support::dtype>(0, "temp");

  auto coords = mesh.coords();

  auto temp_h = Omega_h::HostRead<support::dtype>(temp);
  auto coords_h = Omega_h::HostRead<double>(coords);

  double err2 = 0.0;

  const int nv = mesh.nverts();

  for (int v = 0; v < nv; ++v) {

    double x = coords_h[2 * v]; // x-coordinate

    double exact = 270.0 + 30.0 * x;

    double diff = std::abs(temp_h[v] - exact);

    err2 += diff * diff;
  }

  return std::sqrt(err2 / nv);
}
void PrintTempStats(const Omega_h::Mesh& mesh, const std::string& name, int itr)
{
  auto temp = mesh.get_array<support::dtype>(0, "temp");
  auto temp_h = Omega_h::HostRead<support::dtype>(temp);

  double min_T = std::numeric_limits<double>::max();
  double max_T = -std::numeric_limits<double>::max();
  double sum_T = 0.0;

  for (int i = 0; i < temp_h.size(); ++i) {
    const double T = temp_h[i];
    min_T = std::min(min_T, T);
    max_T = std::max(max_T, T);
    sum_T += T;
  }

  std::cout << "[TEMP_STATS] itr=" << itr << " " << name
            << " size=" << temp_h.size() << " min=" << min_T << " max=" << max_T
            << " mean=" << sum_T / temp_h.size() << "\n";
}
double ComputeAbsoluteError(const mfem::ParMesh& pmesh,
                            const mfem::ParGridFunction& x)
{
  double err2 = 0.0;

  const int nv = pmesh.GetNV();

  for (int v = 0; v < nv; ++v) {

    const double* coord = pmesh.GetVertex(v);

    double xpos = coord[0];

    double exact = 270.0 + 30.0 * xpos;

    double diff = x(v) - exact;

    err2 += diff * diff;
  }

  return std::sqrt(err2 / nv);
}
void initializeFieldWithGids(pcms::FieldT<pcms::Real>* field,
                             pcms::Real multiplier = 1.0)
{
  auto& layout = field->GetLayout();
  auto gids = layout.GetGids();
  const auto n = layout.GetNumOwnedDofHolder();

  Omega_h::HostWrite<pcms::Real> ids(n);
  PCMS_ALWAYS_ASSERT(n == gids.size());
  Kokkos::parallel_for(
    "init gid",
    Kokkos::RangePolicy<pcms::HostMemorySpace::execution_space>(0, n),
    [=](int i) { ids[i] = gids[i] * multiplier; });

  field->SetDOFHolderData(pcms::make_const_array_view(ids));
}
void PrintTempStats(const mfem::ParMesh& pmesh, const mfem::ParGridFunction& x,
                    const std::string& name, int itr)
{
  auto* pfes = x.ParFESpace();

  mfem::Array<int> vdofs;

  double minT = std::numeric_limits<double>::max();
  double maxT = -std::numeric_limits<double>::max();
  double sumT = 0.0;

  const int nv = pmesh.GetNV();

  for (int v = 0; v < nv; ++v) {
    pfes->GetVertexDofs(v, vdofs);

    if (vdofs.Size() == 0)
      continue;

    const int dof = vdofs[0];
    const double T = x(dof);

    minT = std::min(minT, T);
    maxT = std::max(maxT, T);
    sumT += T;
  }

  std::cout << "[TEMP_STATS] itr=" << itr << " " << name << " min=" << minT
            << " max=" << maxT << " mean=" << sumT / nv << "\n";
}
} // namespace support