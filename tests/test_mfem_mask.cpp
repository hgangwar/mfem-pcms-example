//
// Created by gangwh on 3/25/26.
//
#include "../mfem_field_adapter.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>
#include <vector>
#include <sstream>
#include <string>
#include "mfem.hpp"

using pcms::make_array_view;
using pcms::make_const_array_view;
template <typename T>
std::string diagnose_order_or_value_mismatch(const std::vector<T> &expected,
                                             const std::vector<T> &actual,
                                             double tol = 1e-16)
{
  std::ostringstream os;

  if (expected.size() != actual.size()) {
    os << "Size mismatch: expected " << expected.size()
       << ", actual " << actual.size();
    return os.str();
  }

  auto equal_val = [tol](const T &a, const T &b) -> bool {
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

  std::vector<char> used(actual.size(), 0);
  bool same_multiset = true;

  for (int i = 0; i < (int)expected.size(); ++i) {
    bool found = false;
    for (int j = 0; j < (int)actual.size(); ++j) {
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
    for (int i = 0; i < (int)expected.size() && shown < 10; ++i) {
      if (!equal_val(expected[i], actual[i])) {
        os << "  idx " << i
           << ": expected=" << expected[i]
           << ", actual=" << actual[i] << "\n";
        ++shown;
      }
    }
    return os.str();
  }

  os << "Vectors differ in values, not just order.\n";

  int shown = 0;
  for (int i = 0; i < (int)expected.size() && shown < 10; ++i) {
    if (!equal_val(expected[i], actual[i])) {
      os << "  idx " << i
         << ": expected=" << expected[i]
         << ", actual=" << actual[i] << "\n";
      ++shown;
    }
  }

  for (int i = 0; i < (int)expected.size() && shown < 20; ++i) {
    bool found = false;
    for (int j = 0; j < (int)actual.size(); ++j) {
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
  }
  return std::fabs(val1 - val2) < 1e-16;
}

#define TEST_ASSERT(cond, msg)                                    \
  do {                                                            \
    if (!(cond)) {                                                \
      std::cerr << "[FAIL] " << msg << std::endl;                 \
      PCMS_ALWAYS_ASSERT(cond);                                   \
    }                                                             \
  } while (0)

mfem::Vector make_true_gf_data(const mfem::ParGridFunction &gf_data,
                               const mfem::ParFiniteElementSpace &pfes)
{
  mfem::Vector true_gf_data(pfes.GetTrueVSize());

  TEST_ASSERT(gf_data.Size() == true_gf_data.Size(),
              "This test assumes ldofs == tdofs");

  true_gf_data = gf_data;
  return true_gf_data;
}

void set_true_gf_data(mfem::ParGridFunction &gf_data,
                      const mfem::ParFiniteElementSpace &pfes,
                      const mfem::Vector &true_gf_data)
{
  TEST_ASSERT(gf_data.Size() == true_gf_data.Size(),
              "This test assumes ldofs == tdofs");

  gf_data = true_gf_data;
}

void mark_box_attributes(mfem::Mesh &mesh)
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
std::vector<char> build_gid_membership_mask(const std::vector<T> &gids, int size)
{
  std::vector<char> mask(size, 0);
  for (auto gid : gids) {
    TEST_ASSERT(gid >= 0 && gid < size,
                "Adapter returned gid out of bounds for this test");
    mask[static_cast<int>(gid)] = 1;
  }
  return mask;
}

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  {
    const int nx = 30;
    const int ny = 30;
    const double lx = 0.6;
    const double ly = 1.0;
    const int dim = 2;
    const int order = 1;
    const bool use_mask = true;
    const pcms::LO masked_attr = 1;

    mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
        nx, ny, mfem::Element::QUADRILATERAL, true, lx, ly);

    mark_box_attributes(mesh);

    mfem::ParMesh pmesh(MPI_COMM_WORLD, mesh);
    mesh.Clear();

    mfem::H1_FECollection fec(order, dim);
    mfem::ParFiniteElementSpace pfes(&pmesh, &fec);
    mfem::ParGridFunction gf_data(&pfes);

    TEST_ASSERT(pfes.GetTrueVSize() == pmesh.GetNV(),
                "This test assumes H1 order 1 with vertex-aligned true dofs");

    // In this test configuration, ldofs == tdofs, so R and P are identity.
    mfem::Vector true_before(pfes.GetTrueVSize());
    for (int i = 0; i < true_before.Size(); ++i) {
      true_before[i] = 1000.0 + i;
    }
    set_true_gf_data(gf_data, pfes, true_before);

    pcms::MFEMFieldAdapter adapter(
        std::string("mfem_field_adapter"),
        pmesh,
        pfes,
        gf_data,
        use_mask,
        masked_attr);

    std::vector<int> permutation;

    // ============================================================
    // Serialization block
    // ============================================================
    std::cout << "[INFO] Starting serialization test..." << std::endl;

    auto gids_adapter = adapter.GetGids();
    TEST_ASSERT(!gids_adapter.empty(),
                "Serialization failed: adapter returned empty gids");

    std::vector<double> buffer(gids_adapter.size());

    auto packed_size = adapter.Serialize(make_array_view(buffer),
                                         make_const_array_view(permutation));

    TEST_ASSERT((size_t)packed_size == gids_adapter.size(),
                "Serialization failed: packed size does not match adapter gid count");

    auto true_data_before = make_true_gf_data(gf_data, pfes);

    for (int i = 0; i < (int)buffer.size(); ++i) {
      const int gid = static_cast<int>(gids_adapter[i]);
      TEST_ASSERT(gid >= 0 && gid < true_data_before.Size(),
                  "Serialization failed: adapter gid out of bounds");
      TEST_ASSERT(is_close(buffer[i], true_data_before[gid]),
                  "Serialization failed: packed buffer value does not match true data at adapter gid");
    }

    std::cout << "[PASS] Serialization test passed." << std::endl;

    // ============================================================
    // Deserialization block
    // ============================================================
    std::cout << "[INFO] Starting deserialization test..." << std::endl;

    std::vector<double> modified = buffer;
    for (int i = 0; i < (int)modified.size(); ++i) {
      modified[i] = -10.0 * modified[i];
    }

    adapter.Deserialize(make_const_array_view(modified),
                        make_const_array_view(permutation));

    auto true_data_after = make_true_gf_data(gf_data, pfes);
    auto is_masked = build_gid_membership_mask(gids_adapter, true_data_after.Size());

    for (int i = 0; i < true_data_after.Size(); ++i) {
      if (!is_masked[i]) {
        TEST_ASSERT(is_close(true_data_after[i], true_before[i]),
                    "Deserialization failed: unmasked entry was modified");
      }
    }

    std::vector<double> roundtrip(gids_adapter.size());
    auto roundtrip_size = adapter.Serialize(make_array_view(roundtrip),
                                            make_const_array_view(permutation));

    TEST_ASSERT((size_t)roundtrip_size == gids_adapter.size(),
                "Deserialization failed: roundtrip serialize size mismatch");

    std::string diag = diagnose_order_or_value_mismatch(modified, roundtrip);
    TEST_ASSERT(diag == "Vectors match: same values in same order.",
                "Deserialization failed: " << diag);

    for (int i = 0; i < (int)gids_adapter.size(); ++i) {
      const int gid = static_cast<int>(gids_adapter[i]);
      TEST_ASSERT(is_close(true_data_after[gid], roundtrip[i]),
                  "Deserialization failed: true data at adapter gid does not match roundtrip packed value");
    }

    std::cout << "[PASS] Deserialization test passed." << std::endl;
  }

  MPI_Finalize();
  return 0;
}