#include "mfem_adapter_layout.h"
#include "mfem_field_adapter2.h"
#include <vector>
namespace pcms
{

MFEMFieldsAdapterLayout::MFEMFieldsAdapterLayout(
  mfem::ParMesh& pmesh, mfem::ParFiniteElementSpace& pfes,
  mfem::ParGridFunction& gf, std::array<int, 4> nodes_per_dim,
  int num_components, CoordinateSystem coordinate_system, bool use_mask,
  LO attr)
  : pmesh_(pmesh),
    pfes_(pfes),
    gf_(gf),
    num_components_(num_components),
    coordinate_system_(coordinate_system),
    nodes_per_dim_(nodes_per_dim)
{
  PCMS_FUNCTION_TIMER;

  AssertVertexOnly();

  if (use_mask) {
    PCMS_ALWAYS_ASSERT(attr >= 0);
    BuildMaskByAttribute(attr);
  } else {
    has_mask_ = false;
    packed_size_ = pmesh_.GetNV();
  }

  BuildGids();
  BuildCoordinates();
  BuildOwned();
}

void MFEMFieldsAdapterLayout::AssertVertexOnly() const
{
  PCMS_ALWAYS_ASSERT(nodes_per_dim_[0] == 1);
  PCMS_ALWAYS_ASSERT(nodes_per_dim_[1] == 0);
  PCMS_ALWAYS_ASSERT(nodes_per_dim_[2] == 0);
  PCMS_ALWAYS_ASSERT(nodes_per_dim_[3] == 0);

  PCMS_ALWAYS_ASSERT(num_components_ == 1);
  PCMS_ALWAYS_ASSERT(pfes_.GetVDim() == 1);
}

void MFEMFieldsAdapterLayout::BuildMaskByAttribute(LO attr)
{
  const int nverts = pmesh_.GetNV();
  const int ne = pmesh_.GetNE();

  mask_storage_.SetSize(nverts);
  mask_storage_ = 0;

  mfem::Array<int> vert_ids;
  LO count = 0;

  for (int e = 0; e < ne; ++e) {
    if (pmesh_.GetAttribute(e) != attr) {
      continue;
    }

    pmesh_.GetElementVertices(e, vert_ids);

    for (int j = 0; j < vert_ids.Size(); ++j) {
      const int v = vert_ids[j];

      if (mask_storage_[v] == 0) {
        mask_storage_[v] = ++count;
      }
    }
  }

  has_mask_ = true;
  packed_size_ = count;

  mask_view_ = Rank1View<LO, HostMemorySpace>(mask_storage_.GetData(),
                                              mask_storage_.Size());
}

void MFEMFieldsAdapterLayout::BuildGids()
{
  mfem::Array<HYPRE_BigInt> vertex_gids;
  pmesh_.GetGlobalVertexIndices(vertex_gids);

  const int nverts = pmesh_.GetNV();
  const LO nout = has_mask_ ? packed_size_ : nverts;

  gids_host_ = Kokkos::View<GO*, HostMemorySpace>("mfem_gids", nout);

  for (int v = 0; v < nverts; ++v) {
    if (has_mask_ && mask_storage_[v] == 0) {
      continue;
    }

    const LO out = has_mask_ ? mask_storage_[v] - 1 : v;
    gids_host_(out) = static_cast<GO>(vertex_gids[v]);
  }
}

void MFEMFieldsAdapterLayout::BuildCoordinates()
{
  const int nverts = pmesh_.GetNV();
  const int dim = pmesh_.Dimension();
  const LO nout = has_mask_ ? packed_size_ : nverts;

  coords_host_ =
    Kokkos::View<Real**, HostMemorySpace>("mfem_dof_holder_coords", nout, dim);

  mfem::Vector vcoords;
  pmesh_.GetVertices(vcoords);

  for (int i = 0; i < vcoords.Size(); i += dim) {
    const LO v = i / dim;

    if (has_mask_ && mask_storage_[v] == 0) {
      continue;
    }

    const LO out = has_mask_ ? mask_storage_[v] - 1 : v;

    for (int d = 0; d < dim; ++d) {
      coords_host_(out, d) = static_cast<Real>(vcoords[i + d]);
    }
  }
}

void MFEMFieldsAdapterLayout::BuildOwned()
{
  const LO nout = has_mask_ ? packed_size_ : pmesh_.GetNV();

  owned_host_ = Kokkos::View<bool*, HostMemorySpace>("mfem_owned", nout);

  for (LO i = 0; i < nout; ++i) {
    owned_host_(i) = true;
  }
}

std::unique_ptr<FieldT<Real>> MFEMFieldsAdapterLayout::CreateFieldReal() const
{
  return std::make_unique<MFEMFieldAdapter2<Real>>(*this);
}

int MFEMFieldsAdapterLayout::GetNumComponents() const
{
  return num_components_;
}

LO MFEMFieldsAdapterLayout::GetNumOwnedDofHolder() const
{
  return has_mask_ ? packed_size_ : static_cast<LO>(pmesh_.GetNV());
}

GO MFEMFieldsAdapterLayout::GetNumGlobalDofHolder() const
{
  return static_cast<GO>(GetNumOwnedDofHolder());
}

Rank1View<const bool, HostMemorySpace> MFEMFieldsAdapterLayout::GetOwned() const
{
  return Rank1View<const bool, HostMemorySpace>(owned_host_.data(),
                                                owned_host_.extent(0));
}

GlobalIDView<HostMemorySpace> MFEMFieldsAdapterLayout::GetGids() const
{
  return GlobalIDView<HostMemorySpace>(gids_host_.data(), gids_host_.extent(0));
}
CoordinateView<HostMemorySpace>
MFEMFieldsAdapterLayout::GetDOFHolderCoordinates() const
{
  Rank2View<const Real, HostMemorySpace> coords_view(
    coords_host_.data(), coords_host_.extent(0), coords_host_.extent(1));

  return CoordinateView<HostMemorySpace>{coordinate_system_, coords_view};
}

bool MFEMFieldsAdapterLayout::IsDistributed()
{
  return true;
}

EntOffsetsArray MFEMFieldsAdapterLayout::GetEntOffsets() const
{
  EntOffsetsArray offsets{};

  LO offset = 0;

  for (size_t dim = 0; dim < offsets.size(); ++dim) {
    offsets[dim] = offset;

    if (dim <= static_cast<size_t>(pmesh_.Dimension()) && nodes_per_dim_[dim]) {

      if (dim == 0) {
        offset += GetNumOwnedDofHolder();
      } else {
        PCMS_ALWAYS_ASSERT(
          false &&
          "MFEMFieldsAdapterLayout currently supports vertex-only fields.");
      }
    }
  }

  offsets[offsets.size() - 1] = offset;

  return offsets;
}

ReversePartitionMap2 MFEMFieldsAdapterLayout::GetReversePartitionMap(
  const Partition& partition) const
{
  PCMS_FUNCTION_TIMER;

  ReversePartitionMap2 reverse_partition;

  mfem::Vector vcoords;
  const LO dim = pmesh_.Dimension();

  pmesh_.GetVertices(vcoords);

  int local_index = 0;
  std::array<double, 3> coord{0.0, 0.0, 0.0};

  for (int i = 0; i < vcoords.Size(); i += dim) {
    const LO vertex_id = i / dim;

    if (!HasMask() || mask_view_(vertex_id) > 0) {
      std::copy(vcoords.begin() + i, vcoords.begin() + i + dim, coord.begin());

      const auto dr = partition.GetDr(local_index, dim, coord);

      reverse_partition[dr].indices.emplace_back(local_index);

      ++local_index;
    }
  }

  int counter = 0;
  for (const auto& [rank, mapping] : reverse_partition) {
    std::printf("Vertex in Reverse Partition map, rank %d : %zu\n", rank,
                mapping.indices.size());
    counter += static_cast<int>(mapping.indices.size());
  }

  PCMS_ALWAYS_ASSERT(counter == GetNumOwnedDofHolder());

  return reverse_partition;
}

std::array<int, 4> MFEMFieldsAdapterLayout::GetNodesPerDim() const
{
  return nodes_per_dim_;
}

mfem::ParMesh& MFEMFieldsAdapterLayout::GetMesh() const
{
  return pmesh_;
}

mfem::ParFiniteElementSpace& MFEMFieldsAdapterLayout::GetFESpace() const
{
  return pfes_;
}

mfem::ParGridFunction& MFEMFieldsAdapterLayout::GetGridFunction() const
{
  return gf_;
}

bool MFEMFieldsAdapterLayout::HasMask() const
{
  return has_mask_;
}

LO MFEMFieldsAdapterLayout::PackedSize() const
{
  return packed_size_;
}

Rank1View<const LO, HostMemorySpace> MFEMFieldsAdapterLayout::GetMask() const
{
  PCMS_ALWAYS_ASSERT(has_mask_);

  return Rank1View<const LO, HostMemorySpace>(mask_storage_.GetData(),
                                              mask_storage_.Size());
}

} // namespace pcms