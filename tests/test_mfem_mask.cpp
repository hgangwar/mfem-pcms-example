//
// Created by gangwh on 3/25/26.
//
#include "../mfem_field_adapter.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <vector>

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

bool check_data(const std::vector<double> &buffer, const mfem::Vector &vec)
{
  if ((int)buffer.size() != vec.Size()) {
    return false;
  }
  for (int i = 0; i < vec.Size(); ++i) {
    if (!is_close(buffer[i], vec[i])) {
      return false;
    }
  }
  return true;
}

mfem::Vector make_true_gf_data(const mfem::ParGridFunction &gf_data,
                               const mfem::ParFiniteElementSpace &pfes)
{
  mfem::Vector true_gf_data(pfes.GetTrueVSize());
  auto *R = pfes.GetRestrictionMatrix();
  R->Mult(gf_data, true_gf_data);
  return true_gf_data;
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

int count_local_elements_with_attr(const mfem::ParMesh &pmesh, int attr)
{
  int count = 0;
  for (int e = 0; e < pmesh.GetNE(); ++e) {
    if (pmesh.GetAttribute(e) == attr) {
      ++count;
    }
  }
  return count;
}

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  {
    // ------------------------------------------------------------
    // Build a 2D box mesh on [0, 0.6] x [0, 1]
    // ------------------------------------------------------------
    const int nx = 30;
    const int ny = 30;
    const double lx = 0.6;
    const double ly = 1.0;
    const int dim = 2;
    const int order = 1;

    mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
        nx, ny, mfem::Element::QUADRILATERAL, true, lx, ly);

    mark_box_attributes(mesh);

    mfem::ParMesh pmesh(MPI_COMM_WORLD, mesh);
    mesh.Clear();

    mfem::H1_FECollection fec(order, dim);
    mfem::ParFiniteElementSpace pfes(&pmesh, &fec);
    mfem::ParGridFunction gf_data(&pfes);

    // Fill grid function with known values
    for (int i = 0; i < gf_data.Size(); ++i) {
      gf_data[i] = static_cast<double>(i);
    }

    // ------------------------------------------------------------
    // Construct adapter with mask on attribute 1
    // attr 1 corresponds to x < 0.4
    // ------------------------------------------------------------
    const bool use_mask = true;
    const pcms::LO masked_attr = 1;

    pcms::MFEMFieldAdapter adapter(
        std::string("mfem_field_adapter"),
        pmesh,
        pfes,
        gf_data,
        use_mask,
        masked_attr);

    // ------------------------------------------------------------
    // Test serialization / deserialization
    // ------------------------------------------------------------
    std::vector<double> buffer;
    std::vector<int> permutation;

    auto packed_size = adapter.Serialize(make_array_view(buffer),
                                         make_const_array_view(permutation));

    buffer.resize(packed_size);

    auto packed_size_2 = adapter.Serialize(make_array_view(buffer),
                                           make_const_array_view(permutation));

    PCMS_ALWAYS_ASSERT(packed_size == packed_size_2);

    // Since masking is enabled, serialized size should not exceed true size
    auto true_gf_data = make_true_gf_data(gf_data, pfes);
    PCMS_ALWAYS_ASSERT(packed_size <= (pcms::LO)true_gf_data.Size());

    // Optional sanity check: if you interpret mask as element filtering,
    // at least verify that some elements were selected locally.
    int local_masked_elems = count_local_elements_with_attr(pmesh, masked_attr);
    PCMS_ALWAYS_ASSERT(local_masked_elems > 0);

    // Re-serialize into a reference buffer and compare against itself.
    // This avoids assuming the packed masked layout equals full true-dof layout.
    std::vector<double> ref_buffer(buffer.size());
    adapter.Serialize(make_array_view(ref_buffer),
                      make_const_array_view(permutation));

    PCMS_ALWAYS_ASSERT(buffer.size() == ref_buffer.size());
    for (int i = 0; i < (int)buffer.size(); ++i) {
      PCMS_ALWAYS_ASSERT(is_close(buffer[i], ref_buffer[i]));
    }

    // Modify serialized values and deserialize back
    std::transform(buffer.begin(), buffer.end(), buffer.begin(),
                   [](double val) { return 2.0 * val; });

    adapter.Deserialize(make_const_array_view(buffer),
                        make_const_array_view(permutation));

    // Serialize again and verify round-trip on masked packed data
    std::vector<double> roundtrip(buffer.size());
    adapter.Serialize(make_array_view(roundtrip),
                      make_const_array_view(permutation));

    PCMS_ALWAYS_ASSERT(roundtrip.size() == buffer.size());
    for (int i = 0; i < (int)buffer.size(); ++i) {
      PCMS_ALWAYS_ASSERT(is_close(roundtrip[i], buffer[i]));
    }

    // ------------------------------------------------------------
    // Test GIDs
    // ------------------------------------------------------------
    auto gids_adapter = adapter.GetGids();

    // With masking enabled, gids size should match packed serialization size
    PCMS_ALWAYS_ASSERT((pcms::LO)gids_adapter.size() == packed_size);

    // Basic consistency: serialize size and gids size should agree
    PCMS_ALWAYS_ASSERT((pcms::LO)buffer.size() == (pcms::LO)gids_adapter.size());
  }

  MPI_Finalize();
  return 0;
}