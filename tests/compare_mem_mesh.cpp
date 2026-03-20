// compare_mfem_omegah_robust.cpp
#include "mfem.hpp"

#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_vtk.hpp>

#include <vector>
#include <array>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>


using std::array;
using std::vector;

static double sqr(double x) { return x*x; }

struct Stats3 {
  double min =  1e300;
  double max = -1e300;
  double sum = 0.0;
  int n = 0;
  void add(double v){ min=std::min(min,v); max=std::max(max,v); sum+=v; n++; }
  double mean() const { return (n>0)? sum/n : 0.0; }
};

// record layouts for gathering
struct VtxRec {
  std::int64_t gid;
  double x, y;
};

struct ElemRec {
  std::int64_t gid;
  std::int64_t vg0, vg1, vg2; // vertex GIDs of the triangle (sorted)
};
static void Bounds(const vector<array<double,2>> &pts, array<double,2> &mn, array<double,2> &mx)
{
  mn = { 1e300, 1e300 };
  mx = {-1e300,-1e300 };
  for (auto &p : pts)
  {
    mn[0] = std::min(mn[0], p[0]); mn[1] = std::min(mn[1], p[1]);
    mx[0] = std::max(mx[0], p[0]); mx[1] = std::max(mx[1], p[1]);
  }
}

// Brute-force nearest neighbor stats: for each p in A, find min ||p-q|| in B.
// O(N^2) but fine for nx=30,ny=30.
static void NearestNeighborMismatch(const vector<array<double,2>> &A,
                                    const vector<array<double,2>> &B,
                                    double &rms, double &maxd)
{
  double se = 0.0;
  maxd = 0.0;

  for (auto &pa : A)
  {
    double best2 = 1e300;
    for (auto &pb : B)
    {
      double d2 = sqr(pa[0]-pb[0]) + sqr(pa[1]-pb[1]);
      if (d2 < best2) best2 = d2;
    }
    double d = std::sqrt(best2);
    se += d*d;
    maxd = std::max(maxd, d);
  }
  rms = std::sqrt(se / std::max<size_t>(1, A.size()));
}

static double TriArea2D(const array<double,2> &a,
                        const array<double,2> &b,
                        const array<double,2> &c)
{
  // area = 0.5 * |(b-a)x(c-a)|
  double x1 = b[0]-a[0], y1 = b[1]-a[1];
  double x2 = c[0]-a[0], y2 = c[1]-a[1];
  return 0.5 * std::abs(x1*y2 - y1*x2);
}

// Extract MFEM points + triangle connectivity
static void ExtractMFEM(const mfem::Mesh &m,
                        vector<array<double,2>> &pts,
                        vector<array<int,3>> &tris)
{
  pts.resize(m.GetNV());
  for (int i=0;i<m.GetNV();i++){
    const double *v = m.GetVertex(i);
    pts[i] = {v[0], v[1]};
  }

  tris.clear();
  tris.reserve(m.GetNE());
  mfem::Array<int> v;
  for (int e=0;e<m.GetNE();e++){
    m.GetElementVertices(e, v);
    MFEM_VERIFY(v.Size()==3, "MFEM mesh not triangles");
    tris.push_back({v[0], v[1], v[2]});
  }
}

// Extract Omega_h points + triangle connectivity
static void ExtractOmegaH(Omega_h::Mesh &m,
                          vector<array<double,2>> &pts,
                          vector<array<int,3>> &tris)
{
  MFEM_VERIFY(m.dim()==2, "Omega_h mesh not 2D");
  int nv = m.nverts();
  pts.resize(nv);

  auto coords = m.coords();
  auto ch = Omega_h::HostRead<Omega_h::Real>(coords);
  for (int i=0;i<nv;i++){
    pts[i] = { (double)ch[i*2+0], (double)ch[i*2+1] };
  }

  auto ev = m.ask_verts_of(2); // non-const in your version
  auto evh = Omega_h::HostRead<Omega_h::LO>(ev);

  int ne = m.nelems();
  tris.resize(ne);
  for (int e=0;e<ne;e++){
    tris[e] = { (int)evh[e*3+0], (int)evh[e*3+1], (int)evh[e*3+2] };
  }
}

static void EdgeAndAreaStats(const vector<array<double,2>> &pts,
                             const vector<array<int,3>> &tris,
                             Stats3 &edge, Stats3 &area)
{
  for (auto &t : tris)
  {
    const auto &a = pts[t[0]];
    const auto &b = pts[t[1]];
    const auto &c = pts[t[2]];

    double ab = std::hypot(a[0]-b[0], a[1]-b[1]);
    double bc = std::hypot(b[0]-c[0], b[1]-c[1]);
    double ca = std::hypot(c[0]-a[0], c[1]-a[1]);

    edge.add(ab); edge.add(bc); edge.add(ca);
    area.add(TriArea2D(a,b,c));
  }
}
static void CompareGIDsOnRoot(const std::vector<VtxRec>& mfemV,
                             const std::vector<VtxRec>& ohV,
                             const std::vector<ElemRec>& mfemE,
                             const std::vector<ElemRec>& ohE)
{
  // ---- vertices: gid -> coord
  std::unordered_map<std::int64_t, std::array<double,2>> mV, oV;
  mV.reserve(mfemV.size()*2);
  oV.reserve(ohV.size()*2);

  for (auto& r : mfemV) mV[r.gid] = {r.x, r.y};
  for (auto& r : ohV)   oV[r.gid] = {r.x, r.y};

  int64_t commonV = 0, onlyM = 0, onlyO = 0;
  double maxCoordDiff = 0.0;

  for (auto& [gid, xy] : mV){
    auto it = oV.find(gid);
    if (it == oV.end()) { onlyM++; continue; }
    commonV++;
    double dx = xy[0] - it->second[0];
    double dy = xy[1] - it->second[1];
    maxCoordDiff = std::max(maxCoordDiff, std::hypot(dx, dy));
  }
  for (auto& [gid, _] : oV){
    if (mV.find(gid) == mV.end()) onlyO++;
  }

  // ---- elements: gid -> (sorted vertex-gid triple)
  auto key3 = [](std::int64_t a, std::int64_t b, std::int64_t c){
    // cheap hash combiner
    std::uint64_t h = 1469598103934665603ull;
    auto mix = [&](std::uint64_t x){ h ^= x; h *= 1099511628211ull; };
    mix((std::uint64_t)a); mix((std::uint64_t)b); mix((std::uint64_t)c);
    return h;
  };

  std::unordered_map<std::int64_t, std::uint64_t> mE, oE;
  mE.reserve(mfemE.size()*2);
  oE.reserve(ohE.size()*2);

  for (auto& e : mfemE) mE[e.gid] = key3(e.vg0, e.vg1, e.vg2);
  for (auto& e : ohE)   oE[e.gid] = key3(e.vg0, e.vg1, e.vg2);

  int64_t commonE = 0, onlyME = 0, onlyOE = 0, diffConn = 0;
  for (auto& [gid, h] : mE){
    auto it = oE.find(gid);
    if (it == oE.end()) { onlyME++; continue; }
    commonE++;
    if (it->second != h) diffConn++;
  }
  for (auto& [gid, _] : oE){
    if (mE.find(gid) == mE.end()) onlyOE++;
  }

  std::cout << "\n=== GID COMPARISON ===\n";
  std::cout << "Vertices:\n";
  std::cout << "  MFEM unique vtx GIDs:   " << mV.size() << "\n";
  std::cout << "  Omega_h unique vtx GIDs:" << oV.size() << "\n";
  std::cout << "  Common: " << commonV << "  only MFEM: " << onlyM << "  only Omega_h: " << onlyO << "\n";
  std::cout << "  Max coord mismatch on common vertex GIDs: " << maxCoordDiff << "\n";

  std::cout << "Elements:\n";
  std::cout << "  MFEM unique elem GIDs:   " << mE.size() << "\n";
  std::cout << "  Omega_h unique elem GIDs:" << oE.size() << "\n";
  std::cout << "  Common: " << commonE << "  only MFEM: " << onlyME << "  only Omega_h: " << onlyOE << "\n";
  std::cout << "  Connectivity mismatches among common elem GIDs: " << diffConn << "\n";
  std::cout << "======================\n";
}

static void ExtractOmegaH_GIDs(Omega_h::Mesh& m,
                              std::vector<VtxRec>& vtx,
                              std::vector<ElemRec>& elems)
{
  MFEM_VERIFY(m.dim()==2, "Omega_h mesh not 2D");

  // vertex coords + gids
  auto coords = m.coords();
  auto ch = Omega_h::HostRead<Omega_h::Real>(coords);

  auto vglobals = m.globals(0);
  auto vgh = Omega_h::HostRead<Omega_h::GO>(vglobals);

  int nv = m.nverts();
  vtx.resize(nv);
  for (int i=0;i<nv;i++){
    vtx[i] = { (std::int64_t)vgh[i], (double)ch[i*2+0], (double)ch[i*2+1] };
  }

  // element gids + connectivity in vertex GIDs
  auto eglobals = m.globals(2);
  auto egh = Omega_h::HostRead<Omega_h::GO>(eglobals);

  auto ev = m.ask_verts_of(2);
  auto evh = Omega_h::HostRead<Omega_h::LO>(ev);

  int ne = m.nelems();
  elems.resize(ne);
  for (int e=0;e<ne;e++){
    int lv0 = (int)evh[e*3+0];
    int lv1 = (int)evh[e*3+1];
    int lv2 = (int)evh[e*3+2];
    std::int64_t vg0 = (std::int64_t)vgh[lv0];
    std::int64_t vg1 = (std::int64_t)vgh[lv1];
    std::int64_t vg2 = (std::int64_t)vgh[lv2];
    if (vg0 > vg1) std::swap(vg0, vg1);
    if (vg1 > vg2) std::swap(vg1, vg2);
    if (vg0 > vg1) std::swap(vg0, vg1);
    elems[e] = { (std::int64_t)egh[e], vg0, vg1, vg2 };
  }
}
template <class T>
static std::vector<T> GatherToRootBytes(const std::vector<T>& local, MPI_Comm comm, int root=0)
{
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  int local_bytes = int(local.size() * sizeof(T));
  std::vector<int> all_bytes(size, 0);

  MPI_Gather(&local_bytes, 1, MPI_INT,
             all_bytes.data(), 1, MPI_INT,
             root, comm);

  std::vector<int> displs(size, 0);
  int total_bytes = 0;
  if (rank == root) {
    for (int i=0;i<size;i++){
      displs[i] = total_bytes;
      total_bytes += all_bytes[i];
    }
  }

  std::vector<T> global;
  if (rank == root) global.resize(size_t(total_bytes / int(sizeof(T))));

  MPI_Gatherv((void*)local.data(), local_bytes, MPI_BYTE,
              (void*)global.data(), all_bytes.data(), displs.data(), MPI_BYTE,
              root, comm);

  return global;
}

static void ExtractMFEM_GIDs(const mfem::ParMesh& pm,
                            std::vector<VtxRec>& vtx,
                            std::vector<ElemRec>& elems)
{
  // Vertex GIDs
  mfem::Array<HYPRE_BigInt> gv;
  pm.GetGlobalVertexIndices(gv); // size == pm.GetNV()
  vtx.resize(pm.GetNV());
  for (int i=0;i<pm.GetNV();i++){
    const double* p = pm.GetVertex(i);
    vtx[i] = { (std::int64_t)gv[i], p[0], p[1] };
  }

  // Element GIDs
  mfem::Array<HYPRE_BigInt> ge;
  pm.GetGlobalElementIndices(ge); // size == pm.GetNE()

  elems.resize(pm.GetNE());
  mfem::Array<int> v;
  for (int e=0;e<pm.GetNE();e++){
    pm.GetElementVertices(e, v);
    MFEM_VERIFY(v.Size()==3, "MFEM ParMesh not triangles");
    std::int64_t vg0 = (std::int64_t)gv[v[0]];
    std::int64_t vg1 = (std::int64_t)gv[v[1]];
    std::int64_t vg2 = (std::int64_t)gv[v[2]];
    // sort for order-independent compare
    if (vg0 > vg1) std::swap(vg0, vg1);
    if (vg1 > vg2) std::swap(vg1, vg2);
    if (vg0 > vg1) std::swap(vg0, vg1);
    elems[e] = { (std::int64_t)ge[e], vg0, vg1, vg2 };
  }
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  {
    int nx=30, ny=30;
    bool write_vtk = true;

    mfem::OptionsParser args(argc, argv);
    args.AddOption(&nx, "-nx", "--nx", "Elements in x on [0,0.6].");
    args.AddOption(&ny, "-ny", "--ny", "Elements in y on [0,1].");
    args.AddOption(&write_vtk, "-vtk", "--vtk", "-no-vtk", "--no-vtk", "Write VTK outputs.");
    args.Parse();
    if (!args.Good())
    {
      if (mfem::Mpi::WorldRank()==0) args.PrintUsage(std::cout);
      MPI_Finalize();
      return 1;
    }
    if (mfem::Mpi::WorldRank()==0) args.PrintOptions(std::cout);

    // MFEM: triangles on [0,0.6]x[0,1]
    mfem::Mesh mfem_mesh = mfem::Mesh::MakeCartesian2D(nx, ny, mfem::Element::TRIANGLE, true, 0.6, 1.0);
    mfem::ParMesh pm(MPI_COMM_WORLD, mfem_mesh);

    // Omega_h: triangles on [0,0.6]x[0,1] (newer API)
    Omega_h::Library lib(&argc, &argv);
    Omega_h::Mesh oh_mesh = Omega_h::build_box(
      lib.world(),
      OMEGA_H_SIMPLEX,
      0.6, 1, 0,
      nx, ny, 0
    );

    // Extract geometry + connectivity
    vector<array<double,2>> ptsM, ptsO;
    vector<array<int,3>> trisM, trisO;
    ExtractMFEM(mfem_mesh, ptsM, trisM);
    ExtractOmegaH(oh_mesh, ptsO, trisO);

    // Bounding boxes
    array<double,2> mnM, mxM, mnO, mxO;
    Bounds(ptsM, mnM, mxM);
    Bounds(ptsO, mnO, mxO);

    if (mfem::Mpi::WorldRank()==0)
    {
      std::cout << "\n=== OVERALL PICTURE COMPARISON ===\n";
      std::cout << "MFEM:   nv=" << ptsM.size() << " ne=" << trisM.size()
                << "  bbox=[(" << mnM[0] << "," << mnM[1] << "),(" << mxM[0] << "," << mxM[1] << ")]\n";
      std::cout << "OmegaH: nv=" << ptsO.size() << " ne=" << trisO.size()
                << "  bbox=[(" << mnO[0] << "," << mnO[1] << "),(" << mxO[0] << "," << mxO[1] << ")]\n";
    }

    // Point-cloud mismatch (bidirectional)
    double rmsMO, maxMO, rmsOM, maxOM;
    NearestNeighborMismatch(ptsM, ptsO, rmsMO, maxMO); // MFEM -> OmegaH
    NearestNeighborMismatch(ptsO, ptsM, rmsOM, maxOM); // OmegaH -> MFEM

    if (mfem::Mpi::WorldRank()==0)
    {
      std::cout << "Nearest-neighbor mismatch:\n";
      std::cout << "  MFEM -> OmegaH: rms=" << rmsMO << "  max=" << maxMO << "\n";
      std::cout << "  OmegaH -> MFEM: rms=" << rmsOM << "  max=" << maxOM << "\n";
    }

    // Edge length / area stats
    Stats3 eM, aM, eO, aO;
    EdgeAndAreaStats(ptsM, trisM, eM, aM);
    EdgeAndAreaStats(ptsO, trisO, eO, aO);

    // ---- GID compare
    std::vector<VtxRec> mfemV_local, ohV_local;
    std::vector<ElemRec> mfemE_local, ohE_local;

    ExtractMFEM_GIDs(pm, mfemV_local, mfemE_local);
    ExtractOmegaH_GIDs(oh_mesh, ohV_local, ohE_local);

    // gather to root
    auto mfemV_all = GatherToRootBytes(mfemV_local, MPI_COMM_WORLD, 0);
    auto ohV_all   = GatherToRootBytes(ohV_local,   MPI_COMM_WORLD, 0);
    auto mfemE_all = GatherToRootBytes(mfemE_local, MPI_COMM_WORLD, 0);
    auto ohE_all   = GatherToRootBytes(ohE_local,   MPI_COMM_WORLD, 0);

    if (mfem::Mpi::WorldRank() == 0) {
      CompareGIDsOnRoot(mfemV_all, ohV_all, mfemE_all, ohE_all);
    }

    if (mfem::Mpi::WorldRank()==0)
    {
      std::cout << "Edge length stats (min/mean/max):\n";
      std::cout << "  MFEM:   " << eM.min << " / " << eM.mean() << " / " << eM.max << "\n";
      std::cout << "  OmegaH: " << eO.min << " / " << eO.mean() << " / " << eO.max << "\n";
      std::cout << "Triangle area stats (min/mean/max):\n";
      std::cout << "  MFEM:   " << aM.min << " / " << aM.mean() << " / " << aM.max << "\n";
      std::cout << "  OmegaH: " << aO.min << " / " << aO.mean() << " / " << aO.max << "\n";
      std::cout << "=================================\n";
      std::cout << "Interpretation:\n"
                << "  * If bboxes match but NN mismatch is nonzero, the meshes have different vertex layouts.\n"
                << "  * If edge/area stats are similar, they are geometrically comparable even if not identical.\n";
    }

    // Optional: write VTK for direct visual overlay
    if (write_vtk && mfem::Mpi::WorldRank()==0)
    {
      // MFEM VTK
      std::ofstream mfem_out("mfem_mesh.vtk");
      mfem_mesh.PrintVTK(mfem_out);

      // Omega_h VTK (parallel writer still works in serial)
      Omega_h::vtk::write_parallel("omegah_mesh_vtk", &oh_mesh);

      std::cout << "\nWrote:\n  mfem_mesh.vtk\n  omegah_mesh_vtk.pvtu (and pieces)\n";
      std::cout << "Open both in ParaView to visually compare.\n";
    }
  }
  MPI_Finalize();
  return 0;
}