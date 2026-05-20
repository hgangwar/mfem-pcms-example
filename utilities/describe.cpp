#include <mfem.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>

#include <iostream>
#include <map>
#include <limits>
#include <string>

struct Stats
{
  int count = 0;
  double xmin = std::numeric_limits<double>::max();
  double xmax = -std::numeric_limits<double>::max();
  double ymin = std::numeric_limits<double>::max();
  double ymax = -std::numeric_limits<double>::max();
};

void summarize_mfem(const std::string& file)
{
  mfem::Mesh mesh(file.c_str(), 1, 1);

  std::map<int, Stats> stats;

  for (int e = 0; e < mesh.GetNE(); ++e) {
    int attr = mesh.GetAttribute(e);

    mfem::Array<int> verts;
    mesh.GetElementVertices(e, verts);

    for (int i = 0; i < verts.Size(); ++i) {
      auto v = mesh.GetVertex(verts[i]);

      auto& s = stats[attr];
      s.count++;

      s.xmin = std::min(s.xmin, v[0]);
      s.xmax = std::max(s.xmax, v[0]);
      s.ymin = std::min(s.ymin, v[1]);
      s.ymax = std::max(s.ymax, v[1]);
    }
  }

  std::cout << "MFEM mesh attribute summary\n";

  for (auto& [key, val] : stats) {
    std::cout << "attr=" << key << " count=" << val.count << " x=[" << val.xmin
              << "," << val.xmax << "]" << " y=[" << val.ymin << "," << val.ymax
              << "]\n";
  }
}

void summarize_oh(const std::string& file, int argc, char** argv)
{
  Omega_h::Library lib(&argc, &argv);

  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(file, lib.world(), &mesh);

  int dim = mesh.dim();
  auto coords = mesh.coords();

  auto tag = mesh.get_array<int>(dim, "domain");

  auto ev2v = mesh.ask_elem_verts();
  int nv = (dim == 2) ? 3 : 4;

  std::map<int, Stats> stats;

  for (Omega_h::LO e = 0; e < mesh.nelems(); ++e) {
    int key = tag[e];

    for (int j = 0; j < nv; ++j) {
      auto v = ev2v[e * nv + j];

      double x = coords[v * dim + 0];
      double y = coords[v * dim + 1];

      auto& s = stats[key];
      s.count++;

      s.xmin = std::min(s.xmin, x);
      s.xmax = std::max(s.xmax, x);
      s.ymin = std::min(s.ymin, y);
      s.ymax = std::max(s.ymax, y);
    }
  }

  std::cout << "Omega_h tag summary\n";

  for (auto& [key, val] : stats) {
    std::cout << "tag=" << key << " count=" << val.count << " x=[" << val.xmin
              << "," << val.xmax << "]" << " y=[" << val.ymin << "," << val.ymax
              << "]\n";
  }
}

int main(int argc, char** argv)
{
  if (argc < 2) {
    std::cout << "usage: mesh_summary mesh.(msh|mesh|osh)\n";
    return 0;
  }

  std::string file = argv[1];

  if (file.find(".osh") != std::string::npos) {
    summarize_oh(file, argc, argv);
  } else if (file.find(".msh") != std::string::npos ||
             file.find(".mesh") != std::string::npos) {
    summarize_mfem(file);
  } else {
    std::cout << "Unknown mesh type\n";
  }

  return 0;
}