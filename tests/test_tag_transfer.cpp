#include "../schwartz_coupling_support.h"

#include <mfem.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_mesh.hpp>

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <tuple>

using namespace mfem;

constexpr double TOL = 1e-12;

// ------------------------------------------------------------
// Utility: centroid key for matching meshes
// ------------------------------------------------------------
struct Key
{
  long long x, y;
  bool operator<(const Key &o) const
  {
    return std::tie(x, y) < std::tie(o.x, o.y);
  }
};

Key make_key(double x, double y)
{
  return {
      (long long)std::llround(x / TOL),
      (long long)std::llround(y / TOL)};
}

// ------------------------------------------------------------
// Domain selection helper
// domain 0 on [x_start, x_end)
// domain 1 elsewhere
// ------------------------------------------------------------
inline int domain_from_x(double x, double x_start, double x_end)
{
  return (x >= x_start && x < x_end) ? 0 : 1;
}

// ------------------------------------------------------------
// MFEM: mark attributes
// attr 1 <-> domain 0
// attr 2 <-> domain 1
// ------------------------------------------------------------
void mark_mfem(mfem::Mesh &mesh, double x_start, double x_end)
{
  for (int e = 0; e < mesh.GetNE(); ++e)
  {
    Array<int> verts;
    mesh.GetElementVertices(e, verts);

    double xc = 0.0;
    for (int i = 0; i < verts.Size(); ++i)
    {
      xc += mesh.GetVertex(verts[i])[0];
    }
    xc /= verts.Size();

    const int domain = domain_from_x(xc, x_start, x_end);
    const int attr   = domain + 1;

    mesh.SetAttribute(e, attr);
  }

  mesh.SetAttributes();
}

// ------------------------------------------------------------
// Omega_h: mark domain tag
// domain 0 on [x_start, x_end)
// domain 1 elsewhere
// ------------------------------------------------------------
void mark_oh(Omega_h::Mesh &mesh, double x_start, double x_end)
{
  const int dim = mesh.dim();
  auto coords = mesh.coords();
  auto ev2v   = mesh.ask_elem_verts();

  const int nv = (dim == 2) ? 3 : 4;

  Omega_h::Write<Omega_h::I32> dom(mesh.nelems());

  for (Omega_h::LO e = 0; e < mesh.nelems(); ++e)
  {
    double xc = 0.0;
    for (int j = 0; j < nv; ++j)
    {
      const auto v = ev2v[e * nv + j];
      xc += coords[v * dim + 0];
    }
    xc /= nv;

    dom[e] = domain_from_x(xc, x_start, x_end);
  }

  if (mesh.has_tag(dim, "domain"))
  {
    mesh.remove_tag(dim, "domain");
  }

  mesh.add_tag<Omega_h::I32>(dim, "domain", 1,
                             Omega_h::Read<Omega_h::I32>(dom));
}

// ------------------------------------------------------------
// Collect MFEM info
// ------------------------------------------------------------
std::map<Key, int> collect_mfem(const mfem::Mesh &mesh)
{
  std::map<Key, int> map;

  for (int e = 0; e < mesh.GetNE(); ++e)
  {
    Array<int> verts;
    mesh.GetElementVertices(e, verts);

    double xc = 0.0, yc = 0.0;
    for (int i = 0; i < verts.Size(); ++i)
    {
      auto v = mesh.GetVertex(verts[i]);
      xc += v[0];
      yc += v[1];
    }

    xc /= verts.Size();
    yc /= verts.Size();

    map[make_key(xc, yc)] = mesh.GetAttribute(e);
  }

  return map;
}

// ------------------------------------------------------------
// Verify
// ------------------------------------------------------------
bool verify(const mfem::Mesh &mfem_mesh, Omega_h::Mesh &oh_mesh)
{
  auto mfem_map = collect_mfem(mfem_mesh);

  const int dim = oh_mesh.dim();
  auto coords = oh_mesh.coords();
  auto ev2v   = oh_mesh.ask_elem_verts();
  auto dom    = oh_mesh.get_array<Omega_h::I32>(dim, "domain");

  const int nv = (dim == 2) ? 3 : 4;
  std::cout<< "Num of elems: "<< oh_mesh.nelems() << "\n";
  int errors = 0;

  for (Omega_h::LO e = 0; e < oh_mesh.nelems(); ++e)
  {
    double xc = 0.0, yc = 0.0;

    for (int j = 0; j < nv; ++j)
    {
      const auto v = ev2v[e * nv + j];
      xc += coords[v * dim + 0];
      yc += coords[v * dim + 1];
    }

    xc /= nv;
    yc /= nv;

    auto it = mfem_map.find(make_key(xc, yc));
    if (it == mfem_map.end())
    {
      std::cout << "Missing element match\n";
      errors++;
      continue;
    }

    const int mfem_attr = it->second;
    const int oh_dom    = dom[e];

    if (mfem_attr != oh_dom + 1)
    {
      std::cout << "Mismatch at (" << xc << "," << yc << ") "
                << "mfem=" << mfem_attr
                << " oh=" << oh_dom
                << " elem="<< e << "\n";
      errors++;
    }
  }

  if (errors == 0)
  {
    std::cout << "PASS\n";
    return true;
  }
  else
  {
    std::cout << "FAIL: " << errors << " mismatches\n";
    return false;
  }
}

void write_oh_mesh(Omega_h::Mesh &mesh, const std::string &path)
{
  Omega_h::binary::write(path, &mesh);
}

void write_mfem_mesh(const mfem::Mesh &mesh, const std::string &path)
{
  std::ofstream os(path);
  if (!os)
  {
    throw std::runtime_error("Failed to open MFEM output: " + path);
  }
  mesh.Print(os);
}

// ------------------------------------------------------------
// main
// ------------------------------------------------------------
int main(int argc, char *argv[])
{
  if (argc != 5)
  {
    std::cout << "Usage: ./test mesh.mesh mesh.osh x_start x_end\n";
    return 1;
  }

  std::string mfem_mesh_file = argv[1];
  std::string osh_file       = argv[2];
  const double x_start       = std::atof(argv[3]);
  const double x_end         = std::atof(argv[4]);

  if (!(x_start < x_end))
  {
    std::cerr << "Error: require x_start < x_end\n";
    return 2;
  }

  // MFEM reads MFEM mesh file
  mfem::Mesh mfem_mesh(mfem_mesh_file.c_str(), 1, 1);
  mark_mfem(mfem_mesh, x_start, x_end);

  // Omega_h reads .osh
  Omega_h::Library lib(&argc, &argv);
  Omega_h::Mesh oh_mesh(&lib);
  Omega_h::binary::read(osh_file, lib.world(), &oh_mesh);

  mark_oh(oh_mesh, x_start, x_end);

  assert(verify(mfem_mesh, oh_mesh));

  // support::remove_oh_tag(oh_mesh, "domain");
  // support::reset_mfem_attributes(mfem_mesh);

  // write_mfem_mesh(mfem_mesh, "box_A.mesh");
  // write_oh_mesh(oh_mesh, "box_A.osh");

  return 0;
}