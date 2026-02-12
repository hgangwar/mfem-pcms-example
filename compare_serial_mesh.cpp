// compare_box_tri_msh_osh.cpp
//
// Usage (MPI or serial):
//   mpirun -n 1 ./compare_box_tri box_tri
//
// Expects:
//   box_tri.msh   (Gmsh mesh, triangles)
//   box_tri.osh   (Omega_h binary mesh)   OR box_tri.iosh (if that's your extension)
//
// What it reports (on rank 0):
//   - bbox, counts
//   - edge length + tri area stats for each mesh
//   - vertex "GID consistency" by *coordinate-key* (NOT by raw gid==gid)
//   - edge set comparison by *coordinate-key*
//
// Build notes:
//   - Link against MFEM (with MPI) and Omega_h.
//   - In CMake, simplest is target_link_libraries(... mfem Omega_h::omega_h)
//   - CLI compile example (adjust include/lib paths as needed):
//       mpicxx -O2 compare_box_tri_msh_osh.cpp -o compare_box_tri \
//              -I$MFEM_INC -L$MFEM_LIB -lmfem \
//              -I$OMEGAH_INC -L$OMEGAH_LIB -lomega_h
//
#include "mfem.hpp"

#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_file.hpp>   // Omega_h::binary::read (API varies by version)
#include <Omega_h_vtk.hpp>

#include <mpi.h>

#include <vector>
#include <array>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <cstdint>
#include <string>

using std::array;
using std::vector;

static double sqr(double x) { return x * x; }

struct Stats3 {
  double min =  1e300;
  double max = -1e300;
  double sum = 0.0;
  long long n = 0;
  void add(double v){ min=std::min(min,v); max=std::max(max,v); sum+=v; n++; }
  double mean() const { return (n>0)? sum/double(n) : 0.0; }
};

static double TriArea2D(const array<double,2> &a,
                        const array<double,2> &b,
                        const array<double,2> &c)
{
  double x1 = b[0]-a[0], y1 = b[1]-a[1];
  double x2 = c[0]-a[0], y2 = c[1]-a[1];
  return 0.5 * std::abs(x1*y2 - y1*x2);
}

static void Bounds(const vector<array<double,2>> &pts, array<double,2> &mn, array<double,2> &mx)
{
  mn = { 1e300, 1e300 };
  mx = {-1e300,-1e300 };
  for (auto &p : pts) {
    mn[0] = std::min(mn[0], p[0]); mn[1] = std::min(mn[1], p[1]);
    mx[0] = std::max(mx[0], p[0]); mx[1] = std::max(mx[1], p[1]);
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

// -------------------- coordinate keys (robust cross-mesh matching) --------------------
struct XYKey {
  long long ix, iy;
  bool operator==(const XYKey& o) const { return ix==o.ix && iy==o.iy; }
};
struct XYKeyHash {
  size_t operator()(const XYKey& k) const {
    size_t h1 = std::hash<long long>{}(k.ix);
    size_t h2 = std::hash<long long>{}(k.iy);
    return h1 ^ (h2 + 0x9e3779b97f4a7c15ull + (h1<<6) + (h1>>2));
  }
};
static XYKey MakeKey(double x, double y, double tol){
  return XYKey{ (long long)llround(x / tol), (long long)llround(y / tol) };
}

struct EdgeKey {
  XYKey a, b;
  bool operator==(EdgeKey const& o) const {
    return a==o.a && b==o.b;
  }
};
struct EdgeKeyHash {
  size_t operator()(EdgeKey const& e) const {
    size_t h = XYKeyHash{}(e.a);
    size_t g = XYKeyHash{}(e.b);
    return h ^ (g + 0x9e3779b97f4a7c15ull + (h<<6) + (h>>2));
  }
};
static inline void Canonicalize(EdgeKey& ek){
  // order-independent edge
  if (ek.b.ix < ek.a.ix || (ek.b.ix==ek.a.ix && ek.b.iy < ek.a.iy)) std::swap(ek.a, ek.b);
}

// -------------------- MPI gather bytes --------------------
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

// -------------------- records --------------------
struct VtxRec {
  std::int64_t gid;
  double x, y;
};

static void CompareVertexGIDsByCoordOnRoot(const std::vector<VtxRec>& mfemV,
                                          const std::vector<VtxRec>& ohV,
                                          double tol)
{
  std::unordered_map<XYKey, std::int64_t, XYKeyHash> m, o;
  m.reserve(mfemV.size()*2);
  o.reserve(ohV.size()*2);

  for (auto& r : mfemV) m[MakeKey(r.x,r.y,tol)] = r.gid;
  for (auto& r : ohV)   o[MakeKey(r.x,r.y,tol)] = r.gid;

  long long common = 0, onlyM = 0, onlyO = 0, gidDiff = 0;
  for (auto& [k, gidM] : m) {
    auto it = o.find(k);
    if (it == o.end()) { onlyM++; continue; }
    common++;
    if (it->second != gidM) gidDiff++;
  }
  for (auto& [k, _] : o) if (!m.count(k)) onlyO++;

  std::cout << "\n=== VERTEX GID COMPARE (MATCH BY COORD KEY) ===\n";
  std::cout << "Coord tolerance: " << tol << "\n";
  std::cout << "Common coords: " << common
            << "  only MFEM coords: " << onlyM
            << "  only Omega_h coords: " << onlyO << "\n";
  std::cout << "Same coord but different GID: " << gidDiff << "\n";
  std::cout << "================================================\n";
}

static void CompareEdgeSetsByCoordOnRoot(const std::unordered_set<EdgeKey,EdgeKeyHash>& mE,
                                        const std::unordered_set<EdgeKey,EdgeKeyHash>& oE)
{
  long long common=0, onlyM=0, onlyO=0;
  for (auto const& e : mE) {
    if (oE.find(e)!=oE.end()) common++;
    else onlyM++;
  }
  for (auto const& e : oE) if (mE.find(e)==mE.end()) onlyO++;

  std::cout << "\n=== EDGE SET COMPARE (MATCH BY COORD KEY) ===\n";
  std::cout << "MFEM edges:   " << (long long)mE.size() << "\n";
  std::cout << "Omega_h edges:" << (long long)oE.size() << "\n";
  std::cout << "Common: " << common << "  only MFEM: " << onlyM << "  only Omega_h: " << onlyO << "\n";
  std::cout << "================================================\n";
}

// -------------------- Extract MFEM (.msh via MFEM reader) --------------------
static void ExtractMFEM_2D_Tris(const mfem::ParMesh& pm,
                               std::vector<array<double,2>>& pts,
                               std::vector<array<int,3>>& tris,
                               std::vector<VtxRec>& vtx_recs_local,
                               std::unordered_set<EdgeKey,EdgeKeyHash>& edge_set_local,
                               double coord_tol)
{
  // points
  pts.resize(pm.GetNV());
  for (int i=0;i<pm.GetNV();i++){
    const double *v = pm.GetVertex(i);
    pts[i] = {v[0], v[1]};
  }

  // global vertex ids (for reporting only; NOT used to match cross-mesh)
  mfem::Array<HYPRE_BigInt> gv;
  pm.GetGlobalVertexIndices(gv);
  vtx_recs_local.resize(pm.GetNV());
  for (int i=0;i<pm.GetNV();i++){
    vtx_recs_local[i] = { (std::int64_t)gv[i], pts[i][0], pts[i][1] };
  }

  // triangles + edges (edge set by coord-key)
  tris.clear();
  tris.reserve(pm.GetNE());
  edge_set_local.clear();

  mfem::Array<int> v;
  for (int e=0;e<pm.GetNE();e++){
    pm.GetElementVertices(e, v);
    MFEM_VERIFY(v.Size()==3, "MFEM mesh not triangles");
    tris.push_back({v[0], v[1], v[2]});

    // add 3 edges using coordinate keys
    int a=v[0], b=v[1], c=v[2];
    XYKey ka = MakeKey(pts[a][0], pts[a][1], coord_tol);
    XYKey kb = MakeKey(pts[b][0], pts[b][1], coord_tol);
    XYKey kc = MakeKey(pts[c][0], pts[c][1], coord_tol);

    EdgeKey e1{ka,kb}; Canonicalize(e1); edge_set_local.insert(e1);
    EdgeKey e2{kb,kc}; Canonicalize(e2); edge_set_local.insert(e2);
    EdgeKey e3{kc,ka}; Canonicalize(e3); edge_set_local.insert(e3);
  }
}

// -------------------- Extract Omega_h (.osh / .iosh) --------------------
static void ExtractOmegaH_2D_Tris(Omega_h::Mesh& m,
                                 std::vector<array<double,2>>& pts,
                                 std::vector<array<int,3>>& tris,
                                 std::vector<VtxRec>& vtx_recs_local,
                                 std::unordered_set<EdgeKey,EdgeKeyHash>& edge_set_local,
                                 double coord_tol)
{
  MFEM_VERIFY(m.dim()==2, "Omega_h mesh not 2D");
  int nv = m.nverts();
  pts.resize(nv);

  // coords
  auto coords = m.coords();
  auto ch = Omega_h::HostRead<Omega_h::Real>(coords);
  for (int i=0;i<nv;i++){
    pts[i] = { (double)ch[i*2+0], (double)ch[i*2+1] };
  }

  // vertex globals
  auto vglobals = m.globals(0);
  auto vgh = Omega_h::HostRead<Omega_h::GO>(vglobals);
  vtx_recs_local.resize(nv);
  for (int i=0;i<nv;i++){
    vtx_recs_local[i] = { (std::int64_t)vgh[i], pts[i][0], pts[i][1] };
  }

  // triangles
  auto ev = m.ask_verts_of(2);
  auto evh = Omega_h::HostRead<Omega_h::LO>(ev);
  int ne = m.nelems();
  tris.resize(ne);
  edge_set_local.clear();

  for (int e=0;e<ne;e++){
    int a = (int)evh[e*3+0];
    int b = (int)evh[e*3+1];
    int c = (int)evh[e*3+2];
    tris[e] = {a,b,c};

    XYKey ka = MakeKey(pts[a][0], pts[a][1], coord_tol);
    XYKey kb = MakeKey(pts[b][0], pts[b][1], coord_tol);
    XYKey kc = MakeKey(pts[c][0], pts[c][1], coord_tol);

    EdgeKey e1{ka,kb}; Canonicalize(e1); edge_set_local.insert(e1);
    EdgeKey e2{kb,kc}; Canonicalize(e2); edge_set_local.insert(e2);
    EdgeKey e3{kc,ka}; Canonicalize(e3); edge_set_local.insert(e3);
  }
}

// -------------------- Omega_h read helper (handles common API variants) --------------------
static Omega_h::Mesh ReadOmegaHMesh(const std::string& fname, Omega_h::Library& lib)
{
  Omega_h::Mesh m;

  // API variants seen across Omega_h versions:
  //  1) Omega_h::binary::read(fname, comm) -> Mesh
  //  2) Omega_h::binary::read(fname, comm, &mesh)
  //  3) Omega_h::binary::read(fname, &mesh)
  //
  // Try the most common one first; if it doesn't compile in your tree,
  // swap to the commented alternative that matches your install.

  // Variant 1:
  m = Omega_h::binary::read(fname, lib.world()); // common in many installs

  // Variant 2 (uncomment if needed):
  // Omega_h::binary::read(fname, lib.world(), &m);

  // Variant 3 (uncomment if needed):
  // Omega_h::binary::read(fname, &m);

  return m;
}

// -------------------- main --------------------
int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  {
    int rank=0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::string base = "/users/gangwh/src/mfem-pcms-example/mesh/box_tri";
    if (argc >= 2) base = argv[1];

    const std::string msh = base + "./";
    std::string osh = base + ".osh";


    // coordinate-key tolerance:
    // use something much smaller than mesh spacing; 1e-12 is fine for exact box grids.
    const double coord_tol = 1e-12;

    // ----- Read MFEM from .msh (Gmsh) -----
    mfem::Mesh serial_m("/users/gangwh/src/mfem-pcms-example/mesh/box_tri.msh", 1, 1); // 1,1 => generate boundary elements, refine? (standard MFEM ctor)
    mfem::ParMesh pm(MPI_COMM_WORLD, serial_m);

    // ----- Read Omega_h from .osh/.iosh -----
    Omega_h::Library lib(&argc, &argv);

    Omega_h::Mesh oh;
    bool read_ok = false;
    {
      std::ifstream f(osh.c_str(), std::ios::binary);
      if (f.good()) { f.close(); oh = ReadOmegaHMesh(osh, lib); read_ok = true; }
    }
    OMEGA_H_CHECK(read_ok);

    // ----- Extract local data -----
    vector<array<double,2>> ptsM, ptsO;
    vector<array<int,3>> trisM, trisO;
    vector<VtxRec> vtxM_local, vtxO_local;
    std::unordered_set<EdgeKey,EdgeKeyHash> edgesM_local, edgesO_local;

    ExtractMFEM_2D_Tris(pm, ptsM, trisM, vtxM_local, edgesM_local, coord_tol);
    ExtractOmegaH_2D_Tris(oh, ptsO, trisO, vtxO_local, edgesO_local, coord_tol);

    // ----- Basic stats (local) -----
    array<double,2> mnM, mxM, mnO, mxO;
    Bounds(ptsM, mnM, mxM);
    Bounds(ptsO, mnO, mxO);

    Stats3 eM, aM, eO, aO;
    EdgeAndAreaStats(ptsM, trisM, eM, aM);
    EdgeAndAreaStats(ptsO, trisO, eO, aO);

    // ----- Gather vertex records to root -----
    auto vtxM_all = GatherToRootBytes(vtxM_local, MPI_COMM_WORLD, 0);
    auto vtxO_all = GatherToRootBytes(vtxO_local, MPI_COMM_WORLD, 0);

    // ----- Gather edges to root -----
    // pack EdgeKey into a trivially-copyable struct
    struct EdgeRec { long long ax, ay, bx, by; };
    vector<EdgeRec> edgesM_packed, edgesO_packed;
    edgesM_packed.reserve(edgesM_local.size());
    edgesO_packed.reserve(edgesO_local.size());
    for (auto const& e : edgesM_local) edgesM_packed.push_back({e.a.ix,e.a.iy,e.b.ix,e.b.iy});
    for (auto const& e : edgesO_local) edgesO_packed.push_back({e.a.ix,e.a.iy,e.b.ix,e.b.iy});

    auto edgesM_all_p = GatherToRootBytes(edgesM_packed, MPI_COMM_WORLD, 0);
    auto edgesO_all_p = GatherToRootBytes(edgesO_packed, MPI_COMM_WORLD, 0);

    if (rank==0)
    {
      std::cout << "\n=== FILES ===\n";
      std::cout << "MFEM (.msh): " << msh << "\n";

      std::cout << "\n=== BASIC COUNTS (rank-local, but rank0 in serial) ===\n";
      std::cout << "MFEM:   nv=" << ptsM.size() << " ne=" << trisM.size()
                << " bbox=[(" << mnM[0] << "," << mnM[1] << "),(" << mxM[0] << "," << mxM[1] << ")]\n";
      std::cout << "OmegaH: nv=" << ptsO.size() << " ne=" << trisO.size()
                << " bbox=[(" << mnO[0] << "," << mnO[1] << "),(" << mxO[0] << "," << mxO[1] << ")]\n";

      std::cout << "\n=== GEOMETRY STATS ===\n";
      std::cout << "Edge length stats (min/mean/max):\n";
      std::cout << "  MFEM:   " << eM.min << " / " << eM.mean() << " / " << eM.max << "\n";
      std::cout << "  OmegaH: " << eO.min << " / " << eO.mean() << " / " << eO.max << "\n";
      std::cout << "Triangle area stats (min/mean/max):\n";
      std::cout << "  MFEM:   " << aM.min << " / " << aM.mean() << " / " << aM.max << "\n";
      std::cout << "  OmegaH: " << aO.min << " / " << aO.mean() << " / " << aO.max << "\n";

      // ----- Compare vertex GIDs by coordinate-key (the only meaningful cross-file check) -----
      CompareVertexGIDsByCoordOnRoot(vtxM_all, vtxO_all, coord_tol);

      // ----- Rebuild edge sets on root and compare -----
      std::unordered_set<EdgeKey,EdgeKeyHash> edgesM_root, edgesO_root;
      edgesM_root.reserve(edgesM_all_p.size()*2);
      edgesO_root.reserve(edgesO_all_p.size()*2);

      for (auto const& r : edgesM_all_p) {
        EdgeKey ek{ XYKey{r.ax,r.ay}, XYKey{r.bx,r.by} };
        Canonicalize(ek);
        edgesM_root.insert(ek);
      }
      for (auto const& r : edgesO_all_p) {
        EdgeKey ek{ XYKey{r.ax,r.ay}, XYKey{r.bx,r.by} };
        Canonicalize(ek);
        edgesO_root.insert(ek);
      }
      CompareEdgeSetsByCoordOnRoot(edgesM_root, edgesO_root);

      std::cout << "\nInterpretation:\n";
      std::cout << "  - If vertex coords match (Common coords == expected) but 'Same coord but different GID' > 0,\n";
      std::cout << "    then both files have the same geometry but different internal numbering (normal).\n";
      std::cout << "  - If edge sets differ, the triangulation/topology differs (e.g., different diagonal split).\n";
      std::cout << "============================================\n";
    }
  }
  MPI_Finalize();
  return 0;
}
