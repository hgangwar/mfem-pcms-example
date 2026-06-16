#include "../include/mfem_adapter_layout.h"
#include "../include/mfem_field_adapter2.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include "mfem.hpp"

using pcms::make_array_view;
using pcms::make_const_array_view;

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
  for (int i = 0; i < static_cast<int>(expected.size()); ++i) {
    if (!equal_val(expected[i], actual[i])) {
      same_order = false;
      break;
    }
  }

  if (same_order) {
    return "Vectors match: same values in same order.";
  }

  std::vector<char> used(actual.size(), 0);
  bool same_multiset = true;

  for (int i = 0; i < static_cast<int>(expected.size()); ++i) {
    bool found = false;
    for (int j = 0; j < static_cast<int>(actual.size()); ++j) {
      if (!used[j] && equal_val(expected[i], actual[j])) {
        used[j] = 1;
        found = true;
        break;
      }
    }
    if (!found) {
      same_multiset = false;
      break;
    }
  }

  if (same_multiset) {
    os << "Vectors contain the same values, but in different order.\n";
    os << "First few positional mismatches:\n";
    int shown = 0;

    for (int i = 0; i < static_cast<int>(expected.size()) && shown < 10; ++i) {
      if (!equal_val(expected[i], actual[i])) {
        os << "  idx " << i << ": expected=" << expected[i]
           << ", actual=" << actual[i] << "\n";
        ++shown;
      }
    }

    return os.str();
  }

  os << "Vectors differ in values, not just order.\n";

  int shown = 0;
  for (int i = 0; i < static_cast<int>(expected.size()) && shown < 10; ++i) {
    if (!equal_val(expected[i], actual[i])) {
      os << "  idx " << i << ": expected=" << expected[i]
         << ", actual=" << actual[i] << "\n";
      ++shown;
    }
  }

  for (int i = 0; i < static_cast<int>(expected.size()) && shown < 20; ++i) {
    bool found = false;
    for (int j = 0; j < static_cast<int>(actual.size()); ++j) {
      if (equal_val(expected[i], actual[j])) {
        found = true;
        break;
      }
    }

    if (!found) {
      os << "  expected value not found in actual: " << expected[i] << "\n";
      ++shown;
    }
  }

  return os.str();
}

template <typename T, typename T2>
bool is_close(T val1, T2 val2)
{
  if constexpr (std::is_integral_v<T>) {
    return val1 == val2;
  } else {
    return std::fabs(val1 - val2) < 1e-16;
  }
}

#define TEST_ASSERT(cond, msg)                                               \
    do {                                                                     \
    if (!(cond)) {                                                           \
    std::ostringstream _test_assert_os;                                      \
    _test_assert_os << msg;                                                  \
    std::cerr << "[FAIL] " << _test_assert_os.str() << std::endl;            \
    PCMS_ALWAYS_ASSERT(cond);                                                \
    }                                                                        \
} while (0)
mfem::Vector make_true_gf_data(const mfem::ParGridFunction& gf_data,
                               const mfem::ParFiniteElementSpace& pfes)
{
  mfem::Vector true_gf_data(pfes.GetTrueVSize());

  TEST_ASSERT(gf_data.Size() == true_gf_data.Size(),
              "This test assumes ldofs == tdofs");

  true_gf_data = gf_data;
  return true_gf_data;
}

void set_true_gf_data(mfem::ParGridFunction& gf_data,
                      const mfem::ParFiniteElementSpace& pfes,
                      const mfem::Vector& true_gf_data)
{
  TEST_ASSERT(gf_data.Size() == true_gf_data.Size(),
              "This test assumes ldofs == tdofs");

  gf_data = true_gf_data;
}

// Optional: only use this if the input .mesh does not already contain attrs.
// For a box [0,0.6] x [0,1], this marks attr 1 on x < 0.4 and attr 2 otherwise.
void mark_box_attributes(mfem::Mesh& mesh)
{
  for (int e = 0; e < mesh.GetNE(); ++e) {
    mfem::Array<int> verts;
    mesh.GetElementVertices(e, verts);

    double xc = 0.0;
    for (int i = 0; i < verts.Size(); ++i) {
      xc += mesh.GetVertex(verts[i])[0];
    }

    xc /= verts.Size();
    mesh.SetAttribute(e, (xc < 0.4) ? 1 : 2);
  }

  mesh.SetAttributes();
}

template <typename T>
std::vector<char> build_gid_membership_mask(const std::vector<T>& gids,
                                            int size)
{
  std::vector<char> mask(size, 0);

  for (auto gid : gids) {
    TEST_ASSERT(gid >= 0 && gid < size,
                "Adapter returned gid out of bounds for this test");
    mask[static_cast<int>(gid)] = 1;
  }

  return mask;
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

bool parse_bool_arg(const std::string& arg, bool& value)
{
  if (arg == "true" || arg == "1" || arg == "on" || arg == "yes") {
    value = true;
    return true;
  }

  if (arg == "false" || arg == "0" || arg == "off" || arg == "no") {
    value = false;
    return true;
  }

  return false;
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int return_code = 0;

  {
    if (argc < 3) {
      if (rank == 0) {
        std::cerr << "Usage: " << argv[0]
                  << " input.mesh <true|false>\n";
      }

      MPI_Finalize();
      return 1;
    }

    bool use_mask = false;
    const std::string arg = argv[2];

    if (!parse_bool_arg(arg, use_mask)) {
      if (rank == 0) {
        std::cerr << "Invalid boolean argument: " << arg << "\n";
        std::cerr << "Expected one of: true, false, 1, 0, on, off, yes, no\n";
      }

      MPI_Finalize();
      return 1;
    }

    if (rank == 0) {
      std::cout << "[INFO] input mesh = " << argv[1] << "\n";
      std::cout << "[INFO] use_mask = " << std::boolalpha << use_mask << "\n";
    }

    const int dim = 2;
    const int order = 1;
    const pcms::LO masked_attr = 1;

    // ------------------------------------------------------------
    // Read MFEM .mesh
    // ------------------------------------------------------------
    mfem::Mesh mesh(argv[1], 1, 1);

    // Uncomment only if the input .mesh does not already have attrs.
    // mark_box_attributes(mesh);

    mfem::ParMesh pmesh(MPI_COMM_WORLD, mesh);
    mesh.Clear();

    mfem::H1_FECollection fec(order, dim);
    mfem::ParFiniteElementSpace pfes(&pmesh, &fec);
    mfem::ParGridFunction gf_data(&pfes);

    TEST_ASSERT(pfes.GetTrueVSize() == pmesh.GetNV(),
                "This test assumes H1 order 1 with vertex-aligned true dofs");

    mfem::Vector true_before(pfes.GetTrueVSize());

    for (int i = 0; i < true_before.Size(); ++i) {
      true_before[i] = 1000.0 + i;
    }

    set_true_gf_data(gf_data, pfes, true_before);

    pcms::MFEMFieldsAdapterLayout layout(
        pmesh,
        pfes,
        gf_data,
        {1, 0, 0, 0},
        1,
        pcms::CoordinateSystem::Cartesian,
        use_mask,
        masked_attr);

    auto adapter_ptr = layout.CreateFieldReal();

    auto& adapter =
        *adapter_ptr;
    std::vector<int> permutation;

    // ============================================================
    // GID stats block
    // ============================================================
    auto gids_view = layout.GetGids();

    std::vector<pcms::GO> gids_adapter(
        gids_view.data_handle(),
        gids_view.data_handle() + gids_view.size());

    TEST_ASSERT(!gids_adapter.empty(),
                use_mask ? "GetGids returned empty vector for masked mode"
                         : "GetGids returned empty vector for unmasked mode");

    if (!use_mask) {
      TEST_ASSERT(gids_adapter.size() == static_cast<size_t>(pfes.GetTrueVSize()),
                  "Unmasked mode failed: GetGids size should equal true vector size");
    }

    log_gid_stats(gids_adapter, rank);

    // ============================================================
    // Serialization block
    // ============================================================
    std::cout << "[INFO][rank " << rank << "] Starting serialization test with "
              << (use_mask ? "mask" : "no mask") << "...\n";

    std::vector<double> buffer(gids_adapter.size());

    auto packed_size = adapter.Serialize(
    make_array_view(buffer),
    make_const_array_view(permutation));

    printf(" Serialized size: %zu, Get gids size: %zu\n",
           static_cast<size_t>(packed_size),
           gids_adapter.size());

    TEST_ASSERT(static_cast<size_t>(packed_size) == gids_adapter.size(),
                "Serialization failed: packed size does not match adapter gid count");

    auto true_data_before = make_true_gf_data(gf_data, pfes);

    for (int i = 0; i < static_cast<int>(buffer.size()); ++i) {
      const int gid = static_cast<int>(gids_adapter[i]);

      TEST_ASSERT(gid >= 0 && gid < true_data_before.Size(),
                  "Serialization failed: adapter gid out of bounds");

      TEST_ASSERT(is_close(buffer[i], true_data_before[gid]),
                  "Serialization failed: packed buffer value does not match true data");
    }

    std::cout << "[PASS][rank " << rank << "] Serialization test passed.\n";

    // ============================================================
    // Deserialization block
    // ============================================================
    std::cout << "[INFO][rank " << rank << "] Starting deserialization test with "
              << (use_mask ? "mask" : "no mask") << "...\n";

    std::vector<double> modified = buffer;

    for (int i = 0; i < static_cast<int>(modified.size()); ++i) {
      modified[i] = -10.0 * modified[i];
    }

    adapter.Deserialize(make_const_array_view(modified),
                        make_const_array_view(permutation));

    auto true_data_after = make_true_gf_data(gf_data, pfes);

    auto is_selected =
      build_gid_membership_mask(gids_adapter, true_data_after.Size());

    if (use_mask) {
      // Masked mode:
      // Only entries returned by GetGids() should be modified.
      for (int i = 0; i < true_data_after.Size(); ++i) {
        if (!is_selected[i]) {
          TEST_ASSERT(is_close(true_data_after[i], true_before[i]),
                      "Masked deserialization failed: unmasked entry was modified");
        }
      }
    } else {
      // Unmasked mode:
      // Every true DOF should be selected and therefore eligible for modification.
      for (int i = 0; i < true_data_after.Size(); ++i) {
        TEST_ASSERT(is_selected[i],
                    "Unmasked deserialization failed: not all true DOFs were selected");
      }
    }

    std::vector<double> roundtrip(gids_adapter.size());

    auto roundtrip_size =
      adapter.Serialize(make_array_view(roundtrip),
                        make_const_array_view(permutation));

    TEST_ASSERT(static_cast<size_t>(roundtrip_size) == gids_adapter.size(),
                "Deserialization failed: roundtrip serialize size mismatch");

    std::string diag = diagnose_order_or_value_mismatch(modified, roundtrip);

    TEST_ASSERT(diag == "Vectors match: same values in same order.",
                "Deserialization failed: " << diag);

    for (int i = 0; i < static_cast<int>(gids_adapter.size()); ++i) {
      const int gid = static_cast<int>(gids_adapter[i]);

      TEST_ASSERT(gid >= 0 && gid < true_data_after.Size(),
                  "Deserialization failed: adapter gid out of bounds");

      TEST_ASSERT(is_close(true_data_after[gid], roundtrip[i]),
                  "Deserialization failed: true data at adapter gid does not "
                  "match roundtrip packed value");
    }

    std::cout << "[PASS][rank " << rank << "] Deserialization test passed.\n";

    if (rank == 0) {
      std::cout << "[PASS] MFEMFieldAdapter "
                << (use_mask ? "masked" : "unmasked")
                << " serialization/deserialization test passed.\n";
    }
  }

  MPI_Finalize();
  return return_code;
}