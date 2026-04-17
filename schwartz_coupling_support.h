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

#include <mpi.h>
#include "mfem.hpp"
#include "mfem_field_adapter.h"

#include <pcms/pcms.h>
#include <pcms/create_field.h>    // Updated Omega_h field
#include <pcms/transfer_field2.h> // Field transfer methods
#include <Omega_h_mesh.hpp>
#include <Omega_h_array.hpp>
#include <Omega_h_mark.hpp>
/*
#include <pcms/omega_h_field.h>
#include <pcms/coupler.h>
*/
#include <Omega_h_build.hpp>

#include <Omega_h_file.hpp>
#include <Omega_h_vtk.hpp>

typedef pcms::Real dtype;
using namespace mfem;
namespace support
{
// -----------------------------
// FEMSystem
// -----------------------------
struct FEMSystem
{
  ParMesh* pmesh = nullptr;
  H1_FECollection* fec = nullptr;
  ParFiniteElementSpace* fes = nullptr;

  ParBilinearForm* a = nullptr;
  ParLinearForm* b = nullptr;
  ParGridFunction* x = nullptr;

  Array<int> ess_tdofs; // essential TRUE dofs for elimination
};
// -----------------------------
// System Parameters
// -----------------------------
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
// Build FE system: -div(k grad T)=0
// -----------------------------
static FEMSystem Init_FEMSystem(ParMesh* pmesh, int order, double kappa_val)
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
static long SolveSystem(FEMSystem& sys, const std::string& solver_type,
                        const std::string& prec_type, double rel_tol,
                        int max_iter)
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
// Derive from MFEM::Coefficient
class ExactTempCoeff : public Coefficient
{
public:
  double Eval(ElementTransformation& T, const IntegrationPoint& ip) override
  {
    Vector x;
    T.Transform(ip, x);
    return 270.0 + 30.0 * x[0];
  }
};

// -----------------------------
// Utilities
// -----------------------------
static double DefaultTolX(const Mesh& mesh)
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
// Shift Mesh along x axis
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
struct Coupling
{
  std::string name;                                  // Coupler name
  std::unique_ptr<pcms::Coupler> cpl;                // Coupler instance
  std::vector<std::string> app_names;                // Attached apps
  std::vector<std::string> field_names;              // Attached fields
  std::map<std::string, pcms::Application*> apps;    // App_name -> pointer
  std::map<std::string, pcms::CoupledField*> fields; // App_field -> pointer
  bool isServer = false;
};
#include <Omega_h_mesh.hpp>
#include <Omega_h_array.hpp>
#include <iostream>

Omega_h::Write<Omega_h::I8> create_mask(Omega_h::Mesh& mesh,
                                        const char* tag_name, int tag_value)
{
  // mask is defined on vertices
  Omega_h::Write<Omega_h::I8> mask(mesh.nents(0), 0);

  const int dim = mesh.dim();

  // tag is defined on elements, since it was written with add_tag(dim, ...)
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

  if (b.size() != n)
    throw std::runtime_error("ComputeRMS: size mismatch");

  if (n == 0)
    return 0.0;

  // Copy to host
  Omega_h::HostRead<dtype> ha(a);
  Omega_h::HostRead<dtype> hb(b);

  double sum_sq = 0.0;

  for (int i = 0; i < n; ++i) {
    double diff = ha[i] - hb[i];
    sum_sq += diff * diff;
  }

  return std::sqrt(sum_sq / static_cast<dtype>(n));
}

//--------------------------------------------------------------
// Init_Coupler<TAdapter>
//--------------------------------------------------------------
template <typename Adapter_type>
Coupling Init_Coupler(MPI_Comm comm, const std::string& name,
                      const std::vector<std::string>& app_names,
                      const std::vector<std::string>& field_names,
                      bool isServer, const redev::Partition ptn,
                      Adapter_type& Adapter)
{
  Coupling cp;
  cp.name = name;
  cp.app_names = app_names;
  cp.field_names = field_names;
  cp.isServer = isServer;

  if (app_names.size() != field_names.size())
    throw std::runtime_error(
      "Mismatch: app_names and field_names must be of the same size.");

  cp.cpl = std::make_unique<pcms::Coupler>(name, comm, isServer, ptn);

  for (size_t i = 0; i < app_names.size(); ++i) {
    auto* app = cp.cpl->AddApplication(app_names[i]);
    cp.fields[app_names[i]] = app->AddField(field_names[i], std::move(Adapter));
    cp.apps[app_names[i]] = app;
  }
  return cp;
}
// Extract sorted-by-y vertex samples on x=xline.
// Extract sorted-by-y vertex samples on x = xline from an Omega_h vertex scalar
// field. Assumes `field_v` is vertex-associated: field_v.size() ==
// mesh.nverts().
static std::vector<std::pair<double, double>> ExtractVertexLineTrace(
  const Omega_h::Mesh& mesh, Omega_h::Read<Omega_h::Real> field_v, double xline,
  double tol)
{
  std::vector<std::pair<double, double>> trace;
  trace.reserve(static_cast<std::size_t>(mesh.nverts()));

  // Vertex coordinates: coords is (nverts * dim) Reals in interleaved layout.
  auto const coords = mesh.coords();
  int const dim = mesh.dim();
  OMEGA_H_CHECK(dim >= 2);

  for (int vi = 0; vi < mesh.nverts(); ++vi) {
    double const x = coords[vi * dim + 0];
    double const y = coords[vi * dim + 1];

    if (std::abs(x - xline) <= tol) {
      // field_v[vi] is the vertex scalar value.
      trace.emplace_back(y, static_cast<double>(field_v[vi]));
    }
  }

  std::sort(trace.begin(), trace.end(),
            [](auto const& a, auto const& b) { return a.first < b.first; });

  // Dedup by y (same logic as your MFEM version)
  std::vector<std::pair<double, double>> uniq;
  uniq.reserve(trace.size());
  for (auto const& p : trace) {
    if (uniq.empty() || std::abs(p.first - uniq.back().first) > 10 * tol)
      uniq.push_back(p);
    else
      uniq.back().second = p.second; // overwrite (last wins)
  }
  return uniq;
}
// Apply BC condition on Omega_h mesh based on vertex coords
// trace: vector of (y_coord, dof_value)
// mesh:  Omega_h mesh with a vertex tag `tag_name` (default "temp")
// x_line: x coordinate of the line
// x_tol:  tolerance to decide if a vertex is on the line
// y_tol:  tolerance used to match vertex y to a trace y sample (binning)
void FillTagOnXLineFromTrace(
  Omega_h::Mesh& mesh, const std::vector<std::pair<double, double>>& trace,
  double x_line, double tol, const char* tag_name = "temp")
{
  if (trace.empty())
    return;
  if (tol <= 0.0)
    throw std::runtime_error("y_tol must be > 0");

  // Must be a vertex tag (dim = 0)
  if (!mesh.has_tag(0, tag_name)) {
    throw std::runtime_error(std::string("Mesh missing vertex tag: ") +
                             tag_name);
  }

  const int dim = mesh.dim();
  if (dim < 1)
    throw std::runtime_error("Mesh dim invalid");

  const int nverts = mesh.nverts();
  if (nverts == 0)
    return;

  // Build a y->value lookup via binning
  // key = round(y / y_tol)
  std::unordered_map<long long, double> ybin_to_val;
  ybin_to_val.reserve(trace.size() * 2);

  auto ykey = [&](double y) -> long long { return llround(y / tol); };

  for (const auto& p : trace) {
    ybin_to_val[ykey(p.first)] = p.second;
  }

  // Read coordinates on host
  Omega_h::HostRead<Omega_h::Real> hcoords(
    mesh.coords()); // size = nverts * dim

  // Read existing tag (host) and create a writable device array
  Omega_h::Reals temp_in = mesh.get_array<Omega_h::Real>(0, tag_name);
  Omega_h::HostRead<Omega_h::Real> htemp_in(temp_in);

  Omega_h::Write<Omega_h::Real> temp_out_w(nverts);
  Omega_h::HostWrite<Omega_h::Real> htemp_out(temp_out_w);

  // Start with old values, then overwrite those on x_line
  for (int v = 0; v < nverts; ++v) {
    htemp_out[v] = htemp_in[v];
  }

  // Overwrite along x = x_line
  for (int v = 0; v < nverts; ++v) {
    const double x = hcoords[v * dim + 0];
    if (std::abs(x - x_line) > tol)
      continue;

    const double y = (dim > 1) ? hcoords[v * dim + 1] : 0.0;
    const auto it = ybin_to_val.find(ykey(y));
    if (it != ybin_to_val.end()) {
      htemp_out[v] = it->second;
    }
    // else: no matching trace sample for this y-bin; leave existing value
  }

  // Push back onto the mesh (in-place update of the tag)
  mesh.set_tag(0, tag_name, Omega_h::Reals(temp_out_w));
}
// Apply a y->value trace to boundary attribute bdr_attr by setting boundary
// vertices (order=1).
static void ApplyBoundaryTraceByAttr(
  ParMesh& pmesh, ParGridFunction& gf, int bdr_attr,
  const std::vector<std::pair<double, double>>& trace, double tol)
{
  MFEM_VERIFY(!trace.empty(),
              "Empty trace passed to ApplyBoundaryTraceByAttr.");

  auto lookup = [&](double y) -> double {
    auto it =
      std::lower_bound(trace.begin(), trace.end(), std::make_pair(y, -1e300),
                       [](auto& a, auto& b) { return a.first < b.first; });
    if (it != trace.end() && std::abs(it->first - y) <= 10 * tol)
      return it->second;
    if (it != trace.begin()) {
      auto it2 = std::prev(it);
      if (std::abs(it2->first - y) <= 10 * tol)
        return it2->second;
    }
    MFEM_ABORT("Trace lookup failed: y-grid mismatch between subdomains.");
    return 0.0;
  };

  Array<int> verts;
  for (int be = 0; be < pmesh.GetNBE(); be++) {
    const Element* bel = pmesh.GetBdrElement(be);
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

static void ApplyBoundaryConstantByAttr(ParMesh& pmesh, ParGridFunction& gf,
                                        int bdr_attr, double value)
{
  Array<int> verts;
  for (int be = 0; be < pmesh.GetNBE(); be++) {
    const Element* bel = pmesh.GetBdrElement(be);
    if (bel->GetAttribute() != bdr_attr) {
      continue;
    }

    pmesh.GetBdrElementVertices(be, verts);
    for (int j = 0; j < verts.Size(); j++) {
      gf(verts[j]) = value;
    }
  }
}

// Diagnostics: boundary attribute min/max over boundary vertices
static void ReportBdrAttrStats(const ParMesh& pmesh, const ParGridFunction& T,
                               int bdr_attr, const char* name)
{
  Array<int> verts;
  double Tmin = 1e300, Tmax = -1e300;
  int cnt = 0;

  for (int be = 0; be < pmesh.GetNBE(); be++) {
    const Element* bel = pmesh.GetBdrElement(be);
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

static void ReportTraceStats(const std::vector<std::pair<double, double>>& tr,
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

// -----------------------------
// To help register ParaView fields
// -----------------------------
struct OutputPack
{
  ParaViewDataCollection pvd;
  ParGridFunction exact;
  ParGridFunction err;

  OutputPack(const std::string& collection, ParMesh& pm,
             ParFiniteElementSpace& fes)
    : pvd(collection.c_str(), &pm), exact(&fes), err(&fes)
  {
    pvd.SetDataFormat(VTKFormat::BINARY);
    pvd.SetHighOrderOutput(true);
  }
};

static std::vector<std::pair<double, double>> RelaxTrace(
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
static void SaveFields(OutputPack& out, const FEMSystem& sys, int it)
{
  ExactTempCoeff exact;

  // update exact
  out.exact.ProjectCoefficient(exact);

  // update error = T - T_exact
  out.err = *sys.x;
  out.err -= out.exact;

  // time-series metadata
  out.pvd.SetCycle(it);
  out.pvd.SetTime((double)it);

  // write
  out.pvd.Save();
}
OMEGA_H_DEVICE Omega_h::I8 isModelEntInOverlap(const int dim, const int id)
{
  // the TOMMS generated geometric model has
  // entity IDs that increase with the distance
  // from the magnetic axis
  if (dim == 2 && (id >= 22 && id <= 34)) {
    return 1;
  } else if (dim == 1 && (id >= 21 && id <= 34)) {
    return 1;
  } else if (dim == 0 && (id >= 21 && id <= 34)) {
    return 1;
  }
  return 0;
}

/**
 * Create the tag 'isOverlap' for each mesh vertex whose value is 1 if the
 * vertex is classified on a model entity in the closure of the geometric model
 * faces forming the overlap region; the value is 0 otherwise.
 */
Omega_h::Read<Omega_h::I8> markOverlapMeshEntities(Omega_h::Mesh& mesh)
{
  // transfer vtx classification to host
  auto classIds = mesh.get_array<Omega_h::ClassId>(0, "class_id");
  auto classDims = mesh.get_array<Omega_h::I8>(0, "class_dim");
  auto isOverlap = Omega_h::Write<Omega_h::I8>(classIds.size(), "isOverlap");
  auto markOverlap = OMEGA_H_LAMBDA(int i)
  {
    isOverlap[i] = isModelEntInOverlap(classDims[i], classIds[i]);
  };
  Omega_h::parallel_for(classIds.size(), markOverlap);
  auto isOverlap_r = Omega_h::read(isOverlap);
  mesh.add_tag(0, "isOverlap", 1, isOverlap_r);
  return isOverlap_r;
}
Omega_h::HostRead<Omega_h::I8> markMeshOverlapRegion(Omega_h::Mesh& mesh)
{
  auto isOverlap = markOverlapMeshEntities(mesh);
  return Omega_h::HostRead(isOverlap);
}
void write_oh_mesh(Omega_h::Mesh& mesh, const std::string& path)
{
  Omega_h::binary::write(path, &mesh);
}
// MFEM read
mfem::Mesh read_mfem_mesh(const std::string& path)
{
  return mfem::Mesh(path.c_str(), 1, 1);
}
void reset_mfem_attributes(mfem::Mesh& mesh, int attr = 1)
{
  for (int e = 0; e < mesh.GetNE(); ++e) {
    mesh.SetAttribute(e, attr);
  }

  mesh.SetAttributes(); // rebuild attribute list
}
void remove_oh_tag(Omega_h::Mesh& mesh, const std::string& name)
{
  int dim = mesh.dim();

  if (mesh.has_tag(dim, name)) {
    mesh.remove_tag(dim, name);
  }
}
} // namespace support

#endif // PCMS_MFEM_COUPLING_SCHWARTZ_COUPLING_SUPPORT_H
