//
// Created by gangwh on 3/25/26.
//
#include "../mfem_field_adapter.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <vector>
#include <iostream>

#include "mfem.hpp"

using pcms::make_array_view;
using pcms::make_const_array_view;

template <typename T, typename T2>
bool is_close(T val1, T2 val2)
{
  if constexpr (std::is_integral_v<T>) {
    return val1 == val2;
  }
  return std::fabs(val1 - val2) < 1e-16;
}

#define TEST_ASSERT(cond, msg)                                      \
  do {                                                              \
    if (!(cond)) {                                                  \
      std::cerr << "[FAIL] " << msg << std::endl;                   \
      PCMS_ALWAYS_ASSERT(cond);                                     \
    }                                                               \
  } while (0)

mfem::Vector make_true_gf_data(const mfem::ParGridFunction &gf_data,
                               const mfem::ParFiniteElementSpace &pfes)
{
  mfem::Vector true_gf_data(pfes.GetTrueVSize());

  auto *R = pfes.GetRestrictionMatrix();
  if (R) {
    R->Mult(gf_data, true_gf_data);
  } else {
    TEST_ASSERT(gf_data.Size() == true_gf_data.Size(),
                "Restriction matrix is null, but gf_data size does not match true size");
    true_gf_data = gf_data;
  }
  return true_gf_data;
}

void set_true_gf_data(mfem::ParGridFunction &gf_data,
                      const mfem::ParFiniteElementSpace &pfes,
                      const mfem::Vector &true_gf_data)
{
  auto *P = pfes.GetProlongationMatrix();
  if (P) {
    P->Mult(true_gf_data, gf_data);
  } else {
    TEST_ASSERT(gf_data.Size() == true_gf_data.Size(),
                "Prolongation matrix is null, but gf_data size does not match true size");
    gf_data = true_gf_data;
  }
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
template <typename T1, typename T2>
void dump_gid_comparison(const std::vector<T1>& expected,
                         const std::vector<T2>& actual)
{
  std::cout << "\n===== GID COMPARISON =====\n";
  std::cout << "Index | Expected | Adapter\n";

  const int n = std::max(expected.size(), actual.size());

  for (int i = 0; i < n; ++i)
  {
    std::cout << i << " | ";

    if (i < (int)expected.size())
      std::cout << expected[i];
    else
      std::cout << "-";

    std::cout << " | ";

    if (i < (int)actual.size())
      std::cout << actual[i];
    else
      std::cout << "-";

    if (i < (int)expected.size() &&
        i < (int)actual.size() &&
        expected[i] != actual[i])
    {
      std::cout << "  <-- mismatch";
    }

    std::cout << "\n";
  }

  std::cout << "===========================\n\n";
}
std::vector<int> get_expected_masked_vertex_gids(const mfem::ParMesh &pmesh,
                                                 int attr)
{
  std::set<int> unique_vertices;
  mfem::Array<int> verts;

  for (int e = 0; e < pmesh.GetNE(); ++e) {
    if (pmesh.GetAttribute(e) != attr) { continue; }

    pmesh.GetElementVertices(e, verts);
    for (int i = 0; i < verts.Size(); ++i) {
      unique_vertices.insert(verts[i]); // serial assumption
    }
  }

  return std::vector<int>(unique_vertices.begin(), unique_vertices.end());
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
                "This test assumes H1 order 1 in serial: true dofs must equal number of vertices");

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

    auto expected_gids = get_expected_masked_vertex_gids(pmesh, masked_attr);
    TEST_ASSERT(!expected_gids.empty(),
                "Expected masked vertex gids is empty");

    std::vector<int> permutation;

    // ============================================================
    // Serialization block
    // ============================================================
    std::cout << "[INFO] Starting serialization test..." << std::endl;

    std::vector<double> buffer(expected_gids.size());
    auto packed_size = adapter.Serialize(make_array_view(buffer),
                                         make_const_array_view(permutation));
    printf("Size of True vertex set: %d, Mask size: %d", (int)expected_gids.size(), packed_size);
    TEST_ASSERT((int)packed_size == (int)expected_gids.size(),
                "Serialization failed: packed size does not match expected masked gid count");

    auto gids_adapter = adapter.GetGids();

    dump_gid_comparison(expected_gids, gids_adapter);
    TEST_ASSERT((int)gids_adapter.size() == (int)expected_gids.size(),
                "Serialization failed: adapter gids size does not match expected gid count");

    for (int i = 0; i < (int)expected_gids.size(); ++i) {
      TEST_ASSERT(gids_adapter[i] == expected_gids[i],
                  "Serialization failed: adapter gids do not match expected gids");
    }

    auto true_data_before = make_true_gf_data(gf_data, pfes);

    for (int i = 0; i < (int)buffer.size(); ++i) {
      const int gid = expected_gids[i];
      TEST_ASSERT(is_close(buffer[i], true_data_before[gid]),
                  "Serialization failed: packed buffer value does not match expected true dof value");
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

    std::vector<char> is_masked(true_data_after.Size(), 0);
    for (int gid : expected_gids) {
      TEST_ASSERT(gid >= 0 && gid < true_data_after.Size(),
                  "Deserialization failed: masked gid out of bounds");
      is_masked[gid] = 1;
    }

    int packed_idx = 0;
    for (int i = 0; i < true_data_after.Size(); ++i) {
      if (is_masked[i]) {
        TEST_ASSERT(is_close(true_data_after[i], modified[packed_idx]),
                    "Deserialization failed: masked entry was not updated correctly");
        ++packed_idx;
      } else {
        TEST_ASSERT(is_close(true_data_after[i], true_before[i]),
                    "Deserialization failed: unmasked entry was modified");
      }
    }

    TEST_ASSERT(packed_idx == (int)modified.size(),
                "Deserialization failed: packed index count mismatch");

    std::vector<double> roundtrip(modified.size());
    auto roundtrip_size = adapter.Serialize(make_array_view(roundtrip),
                                            make_const_array_view(permutation));

    TEST_ASSERT((int)roundtrip_size == (int)modified.size(),
                "Deserialization failed: roundtrip serialize size mismatch");

    for (int i = 0; i < (int)modified.size(); ++i) {
      TEST_ASSERT(is_close(roundtrip[i], modified[i]),
                  "Deserialization failed: roundtrip packed data mismatch");
    }

    std::cout << "[PASS] Deserialization test passed." << std::endl;
  }

  MPI_Finalize();
  return 0;
}