//
// Created by gangwh on 4/7/26.
//
#include <redev.h>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <pcms/utility/assert.h>
#include <pcms.h>
using pcms::make_array_view;
using pcms::make_const_array_view;

template <typename T, typename T2>
bool is_close(T a, T2 b)
{
  if constexpr (std::is_integral_v<T>) {
    return a == b;
  }
  return std::fabs(a - b) < 1e-16;
}

#define TEST_ASSERT(cond, msg)                                                 \
  do {                                                                         \
    if (!(cond)) {                                                             \
      std::cerr << "[FAIL] " << msg << std::endl;                              \
      PCMS_ALWAYS_ASSERT(cond);                                                \
    }                                                                          \
  } while (0)

template <typename T>
std::string diagnose_order_or_value_mismatch(const std::vector<T>& expected,
                                             const std::vector<T>& actual,
                                             double tol = 1e-16)
{
  std::ostringstream os;

  if (expected.size() != actual.size()) {
    os << "Size mismatch: expected " << expected.size() << ", actual "
       << actual.size();
    return os.str();
  }

  auto equal_val = [tol](const T& a, const T& b) -> bool {
    if constexpr (std::is_integral_v<T>) {
      return a == b;
    } else {
      return std::fabs(a - b) < tol;
    }
  };

  bool same_order = true;
  for (int i = 0; i < (int)expected.size(); ++i) {
    if (!equal_val(expected[i], actual[i])) {
      same_order = false;
      break;
    }
  }

  if (same_order) {
    return "Vectors match: same values in same order.";
  }

  os << "Vectors differ.\n";
  int shown = 0;
  for (int i = 0; i < (int)expected.size() && shown < 10; ++i) {
    if (!equal_val(expected[i], actual[i])) {
      os << "  idx " << i << ": expected=" << expected[i]
         << ", actual=" << actual[i] << "\n";
      ++shown;
    }
  }
  return os.str();
}

template <typename T>
void log_gid_stats(const std::vector<T>& gids, int rank)
{
  std::ostringstream os;
  os << "[INFO][rank " << rank << "] GetGids stats\n";
  os << "  length = " << gids.size() << "\n";

  if (gids.empty()) {
    os << "  values = <empty>\n";
    std::cout << os.str() << std::flush;
    return;
  }

  auto [min_it, max_it] = std::minmax_element(gids.begin(), gids.end());
  os << "  min = " << *min_it << "\n";
  os << "  max = " << *max_it << "\n";

  const int preview = std::min<int>(10, gids.size());
  os << "  first " << preview << " gids = ";
  for (int i = 0; i < preview; ++i) {
    os << gids[i];
    if (i + 1 < preview) {
      os << ", ";
    }
  }
  os << "\n";

  std::set<T> uniq(gids.begin(), gids.end());
  os << "  unique count = " << uniq.size() << "\n";

  std::cout << os.str() << std::flush;
}

// -----------------------------------------------------------------------------
// Create vertex mask from an element tag written as add_tag(dim, ...)
// mask[v] = 1 if vertex v belongs to at least one tagged element with tag_value
// -----------------------------------------------------------------------------
Omega_h::Write<Omega_h::I8> create_mask(Omega_h::Mesh& mesh,
                                        const char* tag_name, int tag_value)
{
  Omega_h::Write<Omega_h::I8> mask(mesh.nents(0), 0);

  const int dim = mesh.dim();
  auto tag = mesh.get_array<int>(dim, tag_name);
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
    TEST_ASSERT(false, "Unsupported mesh dimension");
  }

  std::cout << "[INFO] Number of unique vertices filtered: "
            << unique_vertices_filtered << std::endl;

  return mask;
}

// -----------------------------------------------------------------------------
// Add a vertex field tag whose value equals the vertex global id.
// This makes serialize-vs-GetGids checking easy.
// If your mesh does not carry vertex globals, replace this with any known
// field.
// -----------------------------------------------------------------------------
template <typename dtype>
void add_gid_field_tag(Omega_h::Mesh& mesh, const std::string& field_name)
{
  auto vgids = mesh.globals(0); // vertex global ids
  TEST_ASSERT(vgids.size() == mesh.nents(0), "vertex globals size mismatch");

  Omega_h::Write<dtype> field(mesh.nents(0));
  for (Omega_h::LO v = 0; v < mesh.nents(0); ++v) {
    field[v] = static_cast<dtype>(vgids[v]);
  }

  if (mesh.has_tag(0, field_name)) {
    mesh.remove_tag(0, field_name);
  }
  mesh.add_tag<dtype>(0, field_name, 1, Omega_h::Read<dtype>(field));
}

// If you want a snapshot of the current vertex field values:
template <typename dtype>
std::vector<dtype> get_vertex_field(const Omega_h::Mesh& mesh,
                                    const std::string& field_name)
{
  auto field = mesh.get_array<dtype>(0, field_name);
  std::vector<dtype> out(field.size());
  for (int i = 0; i < field.size(); ++i) {
    out[i] = field[i];
  }
  return out;
}

template <typename T>
std::vector<char> build_gid_membership_mask(const std::vector<T>& gids,
                                            int size)
{
  std::vector<char> mask(size, 0);
  for (auto gid : gids) {
    TEST_ASSERT(gid >= 0 && gid < size, "gid out of bounds");
    mask[static_cast<int>(gid)] = 1;
  }
  return mask;
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  {
    if (argc < 2) {
      if (rank == 0) {
        std::cerr << "Usage: " << argv[0] << " input.osh\n";
      }
      MPI_Finalize();
      return 1;
    }

    using dtype = double;
    const std::string field_name = "test_field";
    const char* elem_tag_name = "domain";
    const int masked_tag_value = 0;

    Omega_h::Library lib(&argc, &argv);
    Omega_h::Mesh mesh_A(&lib);
    Omega_h::binary::read(argv[1], lib.world(), &mesh_A);

    // mask from element tag -> vertex mask
    auto is_overlap_A = create_mask(mesh_A, elem_tag_name, masked_tag_value);

    // add a vertex field equal to vertex global id
    add_gid_field_tag<dtype>(mesh_A, field_name);

    // construct adapter
    auto adapter_A =
      pcms::OmegaHFieldAdapter<dtype>(field_name, mesh_A, is_overlap_A);

    // Define Partition
    redev::LOs ranks(1);
    std::iota(ranks.begin(), ranks.end(), 0);
    redev::Reals cuts = {0};
    auto dim = mesh_A.dim();
    auto partition = redev::Partition{redev::RCBPtn{dim, ranks, cuts}};

    // ============================================================
    // GID stats block
    // ============================================================
    auto gids_adapter = adapter_A.GetGids();
    const pcms::ReversePartitionMap reverse_partition =
      adapter_A.GetReversePartitionMap(pcms::Partition{partition});
    auto out_message = pcms::ConstructOutMessage(reverse_partition);
    auto permutation = pcms::ConstructPermutation(reverse_partition);

    TEST_ASSERT(!gids_adapter.empty(),
                "GetGids returned empty vector for masked Omega_h adapter");

    log_gid_stats(gids_adapter, rank);

    // ============================================================
    // Serialization block
    // ============================================================
    std::cout << "[INFO][rank " << rank << "] Starting serialization test...\n";

    std::vector<dtype> buffer(gids_adapter.size());

    auto packed_size = adapter_A.Serialize(make_array_view(buffer),
                                           make_const_array_view(permutation));

    TEST_ASSERT((size_t)packed_size == gids_adapter.size(),
                "Serialize size does not match GetGids size");

    // Since field value == vertex global id, serialized values should match
    // gids
    for (int i = 0; i < (int)buffer.size(); ++i) {
      TEST_ASSERT(is_close(buffer[i], static_cast<dtype>(gids_adapter[i])),
                  "Serialized value does not match adapter gid");
    }

    std::cout << "[PASS][rank " << rank << "] Serialization test passed.\n";

    // ============================================================
    // Deserialization block
    // ============================================================
    std::cout << "[INFO][rank " << rank
              << "] Starting deserialization test...\n";

    auto field_before = get_vertex_field<dtype>(mesh_A, field_name);

    std::vector<dtype> modified = buffer;
    for (int i = 0; i < (int)modified.size(); ++i) {
      modified[i] = -10.0 * modified[i];
    }

    adapter_A.Deserialize(make_const_array_view(modified),
                          make_const_array_view(permutation));

    auto field_after = get_vertex_field<dtype>(mesh_A, field_name);
    auto is_masked =
      build_gid_membership_mask(gids_adapter, field_after.size());

    // Unmasked entries should remain unchanged.
    for (int i = 0; i < (int)field_after.size(); ++i) {
      if (!is_masked[i]) {
        TEST_ASSERT(is_close(field_after[i], field_before[i]),
                    "Unmasked vertex field entry was modified");
      }
    }

    // Round-trip serialize should match modified packed buffer.
    std::vector<dtype> roundtrip(gids_adapter.size());
    auto roundtrip_size = adapter_A.Serialize(
      make_array_view(roundtrip), make_const_array_view(permutation));

    TEST_ASSERT((size_t)roundtrip_size == gids_adapter.size(),
                "Roundtrip serialize size mismatch");

    std::string diag = diagnose_order_or_value_mismatch(modified, roundtrip);
    TEST_ASSERT(diag == "Vectors match: same values in same order.",
                "Roundtrip mismatch: " << diag);

    std::cout << "[PASS][rank " << rank << "] Deserialization test passed.\n";
  }

  MPI_Finalize();
  return 0;
}