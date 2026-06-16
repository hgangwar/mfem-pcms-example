#ifndef PCMS_MFEM_COUPLING_MFEM_ADAPTER_LAYOUT_H
#define PCMS_MFEM_COUPLING_MFEM_ADAPTER_LAYOUT_H

#include <pcms/coordinate_system.h>
#include <pcms/field.h>
#include <pcms/field_layout.h>
#include <pcms/partition.h>
#include <pcms/utility/arrays.h>
#include <pcms/utility/assert.h>
#include <pcms/utility/profile.h>
#include <pcms/utility/types.h>

#include <mfem.hpp>

#include <array>
#include <memory>
#include <string>

namespace pcms
{

template <typename T>
class MFEMFieldAdapter2;

class MFEMFieldsAdapterLayout : public FieldLayout
{
public:
  using value_type = Real;

  MFEMFieldsAdapterLayout(mfem::ParMesh& pmesh,
                          mfem::ParFiniteElementSpace& pfes,
                          mfem::ParGridFunction& gf,
                          std::array<int, 4> nodes_per_dim, int num_components,
                          CoordinateSystem coordinate_system,
                          bool use_mask = false, LO attr = -1);

  std::unique_ptr<FieldT<Real>> CreateFieldReal() const override;

  int GetNumComponents() const override;

  LO GetNumOwnedDofHolder() const override;
  GO GetNumGlobalDofHolder() const override;

  Rank1View<const bool, HostMemorySpace> GetOwned() const override;
  GlobalIDView<HostMemorySpace> GetGids() const override;
  CoordinateView<HostMemorySpace> GetDOFHolderCoordinates() const override;

  bool IsDistributed() override;

  EntOffsetsArray GetEntOffsets() const override;

  ReversePartitionMap2 GetReversePartitionMap(const Partition& partition) const;

  std::array<int, 4> GetNodesPerDim() const;

  mfem::ParMesh& GetMesh() const;
  mfem::ParFiniteElementSpace& GetFESpace() const;
  mfem::ParGridFunction& GetGridFunction() const;

  bool HasMask() const;
  LO PackedSize() const;
  Rank1View<const LO, HostMemorySpace> GetMask() const;

private:
  void AssertVertexOnly() const;
  void BuildMaskByAttribute(LO attr);
  void BuildGids();
  void BuildCoordinates();
  void BuildOwned();

private:
  mfem::ParMesh& pmesh_;
  mfem::ParFiniteElementSpace& pfes_;
  mfem::ParGridFunction& gf_;

  int num_components_;
  CoordinateSystem coordinate_system_;
  std::array<int, 4> nodes_per_dim_;

  bool has_mask_ = false;
  LO packed_size_ = 0;

  mfem::Array<LO> mask_storage_;
  Rank1View<LO, HostMemorySpace> mask_view_;

  Kokkos::View<GO*, HostMemorySpace> gids_host_;
  Kokkos::View<Real**, HostMemorySpace> coords_host_;
  Kokkos::View<bool*, HostMemorySpace> owned_host_;
};

} // namespace pcms

#endif