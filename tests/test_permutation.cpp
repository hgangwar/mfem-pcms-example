//
// Created by gangwh on 6/1/26.
//
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>

#include <pcms.h>
#include <pcms/utility/assert.h>

#include "../include/coupling_support.h"
#include "../include/mfem_field_adapter.h"


#define TEST_ASSERT(cond, msg)                                                 \
  do {                                                                         \
    if (!(cond)) {                                                             \
      std::cerr << "[FAIL] " << msg << std::endl;                              \
      PCMS_ALWAYS_ASSERT(cond);                                                \
    }                                                                          \
  } while (0)

template <typename T>
void PrintGids(const std::string& name, const std::vector<T>& gids, int rank)
{
  if (rank != 0) return;

  std::cout << "\n[INFO] " << name << "\n";
  std::cout << "  size = " << gids.size() << "\n";

  if (gids.empty()) return;

  auto [min_it, max_it] = std::minmax_element(gids.begin(), gids.end());
  std::set<T> uniq(gids.begin(), gids.end());

  std::cout << "  min = " << *min_it << "\n";
  std::cout << "  max = " << *max_it << "\n";
  std::cout << "  unique = " << uniq.size() << "\n";

  std::cout << "  first gids = ";
  for (size_t i = 0; i < std::min<size_t>(20, gids.size()); ++i) {
    std::cout << gids[i] << " ";
  }
  std::cout << "\n";
}

template <typename T>
std::string DiagnosePermutationFailure(const std::vector<T>& mfem_gids,
                                       const std::vector<T>& oh_gids)
{
  std::ostringstream os;

  os << "\nMFEM/Omega_h masked GID permutation failed\n";
  os << "  MFEM size    = " << mfem_gids.size() << "\n";
  os << "  Omega_h size = " << oh_gids.size() << "\n";

  std::multiset<T> mfem_set(mfem_gids.begin(), mfem_gids.end());
  std::multiset<T> oh_set(oh_gids.begin(), oh_gids.end());

  auto oh_copy = oh_set;

  int shown = 0;
  os << "\nGIDs in MFEM but missing in Omega_h:\n";
  for (const auto& gid : mfem_set) {
    auto it = oh_copy.find(gid);
    if (it == oh_copy.end()) {
      os << "  " << gid << "\n";
      if (++shown >= 30) break;
    } else {
      oh_copy.erase(it);
    }
  }

  auto mfem_copy = mfem_set;

  shown = 0;
  os << "\nGIDs in Omega_h but missing in MFEM:\n";
  for (const auto& gid : oh_set) {
    auto it = mfem_copy.find(gid);
    if (it == mfem_copy.end()) {
      os << "  " << gid << "\n";
      if (++shown >= 30) break;
    } else {
      mfem_copy.erase(it);
    }
  }

  os << "\nFirst MFEM gids:\n  ";
  for (size_t i = 0; i < std::min<size_t>(30, mfem_gids.size()); ++i)
    os << mfem_gids[i] << " ";

  os << "\nFirst Omega_h gids:\n  ";
  for (size_t i = 0; i < std::min<size_t>(30, oh_gids.size()); ++i)
    os << oh_gids[i] << " ";

  os << "\n";

  return os.str();
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);

  int rank = 0;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);

  if (argc < 3) {
    if (rank == 0) {
      std::cerr << "Usage:\n"
                << "  " << argv[0]
                << " mfem_mesh.mesh omega_mesh.osh [mfem_attr] [omega_tag]\n\n"
                << "Example:\n"
                << "  " << argv[0] << " A.mesh A.osh 1 0\n";
    }

    MPI_Finalize();
    return 1;
  }

  const std::string mfem_mesh_file = argv[1];
  const std::string omega_mesh_file = argv[2];

  const pcms::LO mfem_attr = argc > 3 ? std::stoi(argv[3]) : 1;
  const pcms::LO omega_tag = argc > 4 ? std::stoi(argv[4]) : 0;

  const std::string field_name = "temp";
  constexpr int order = 1;
  constexpr double kappa = 1.0;

  {
    // -------------------------------------------------------------------------
    // Omega_h adapter
    // -------------------------------------------------------------------------

    Omega_h::Library lib(nullptr, nullptr, comm);
    auto world = lib.world();

    Omega_h::Mesh oh_mesh(&lib);
    Omega_h::binary::read(omega_mesh_file, world, &oh_mesh);

    const auto nverts = oh_mesh.nverts();

    Omega_h::Read<support::dtype> init(nverts, 280.0);
    oh_mesh.add_tag<support::dtype>(Omega_h::VERT, field_name, 1, init);

    auto is_overlap_oh = support::create_mask(oh_mesh, "domain", omega_tag);

    auto oh_adapter =
      pcms::OmegaHFieldAdapter<support::dtype>(field_name, oh_mesh,
                                               is_overlap_oh);

    // -------------------------------------------------------------------------
    // MFEM adapter
    // -------------------------------------------------------------------------

    mfem::Mesh serial_mesh(mfem_mesh_file.c_str(), 1, 1);
    mfem::ParMesh pmesh(comm, serial_mesh);

    support::FEMSystem fem = support::Init_FEMSystem(&pmesh, order, kappa);

    *(fem.x) = 280.0;

    const bool use_mask = true;

    auto mfem_adapter =
      pcms::MFEMFieldAdapter(field_name, *fem.pmesh, *fem.fes,
                                             *fem.x, use_mask, mfem_attr);

    // -------------------------------------------------------------------------
    // Test: masked GetGids from both adapters must have same size, same set,
    // and same local packed ordering.
    // -------------------------------------------------------------------------

    auto mfem_gids = mfem_adapter.GetGids();
    auto oh_gids = oh_adapter.GetGids();

    PrintGids("MFEM masked GetGids", mfem_gids, rank);
    PrintGids("Omega_h masked GetGids", oh_gids, rank);

    TEST_ASSERT(!mfem_gids.empty(), "MFEM masked GetGids returned empty.");
    TEST_ASSERT(!oh_gids.empty(), "Omega_h masked GetGids returned empty.");

    TEST_ASSERT(mfem_gids.size() == oh_gids.size(),
                "Masked GID size mismatch. MFEM size = "
                  << mfem_gids.size()
                  << ", Omega_h size = " << oh_gids.size());

    // First test: same GID set, regardless of order.
    const bool same_gid_set =
      std::is_permutation(mfem_gids.begin(), mfem_gids.end(), oh_gids.begin());

    if (!same_gid_set && rank == 0) {
      std::cerr << DiagnosePermutationFailure(mfem_gids, oh_gids) << std::endl;
    }

    TEST_ASSERT(same_gid_set,
                "MFEM/Omega_h masked GetGids are not permutations.");

    // Second test: same packed local ordering.
    const bool same_local_order = (mfem_gids == oh_gids);

    if (!same_local_order && rank == 0) {
      std::cerr << "\nMFEM/Omega_h masked GID local-order mismatch\n";

      int shown = 0;
      for (size_t i = 0; i < mfem_gids.size(); ++i) {
        if (mfem_gids[i] != oh_gids[i]) {
          std::cerr << "  packed_idx " << i
                    << " : MFEM gid = " << mfem_gids[i]
                    << ", Omega_h gid = " << oh_gids[i]
                    << "\n";

          if (++shown >= 30) break;
        }
      }
    }

    TEST_ASSERT(same_local_order,
                "MFEM/Omega_h masked GetGids have different packed ordering.");

    if (rank == 0) {
      std::cout << "\n[PASS] MFEM/Omega_h masked GetGids are permutations.\n";
      std::cout << "[PASS] MFEM/Omega_h masked GetGids have same packed ordering.\n";
    }
    support::DestroyFEMSystem(fem);
  }

  MPI_Finalize();
  return 0;
}