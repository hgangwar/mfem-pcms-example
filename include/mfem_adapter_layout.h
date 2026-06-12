//
// Created by gangwh on 6/12/26.
//

#ifndef PCMS_MFEM_COUPLING_MFEM_ADAPTER_LAYOUT_H
#define PCMS_MFEM_COUPLING_MFEM_ADAPTER_LAYOUT_H

#include "pcms/field_layout.h"
#include "pcms/coordinate_system.h"
#include "pcms/field.h"

#include <array>
#include "pcms/utility/arrays.h"
#include "mfem.hpp"
namespace pcms
{
class MfemAdapterLayout : public FieldLayout
{
public:
  using value_type = Real;
  MfemAdapterLayout(mfem::ParMesh& pmesh, std::array<int, 4> nodes_per_dim,
                          int num_components,
                          CoordinateSystem coordinate_system,
                          std::string global_id_name = "global");

  std::unique_ptr<FieldT<Real>> CreateFieldReal() const override;

  int GetNumComponents() const override;
  // nodes for standard lagrange FEM
  LO GetNumOwnedDofHolder() const override;
  GO GetNumGlobalDofHolder() const override;

  Rank1View<const bool, HostMemorySpace> GetOwned() const override;
  GlobalIDView<HostMemorySpace> GetGids() const override;
  CoordinateView<HostMemorySpace> GetDOFHolderCoordinates() const override;

  // returns true if the field layout is distributed
  // if the field layout is distributed, the owned and global dofs are the same
  bool IsDistributed() override;

  EntOffsetsArray GetEntOffsets() const override;

  ReversePartitionMap2 GetReversePartitionMap(
    const redev::Partition& partition) const override;

  std::array<int, 4> GetNodesPerDim() const;
  size_t GetNumEnts() const;
  mfem::ParMesh& GetMesh() const;

private:

  mfem::ParMesh& mesh_;
  mfem::Array<pcms::LO> gids_;
  std::string global_id_name_;
  int num_components_;
  CoordinateSystem coordinate_system_;
  std::array<int, 4> nodes_per_dim_;
  Kokkos::View<Real**> dof_holder_coords_;
  Kokkos::View<Real**, HostMemorySpace> dof_holder_coords_host_;
  Kokkos::View<bool*> owned_;
  Kokkos::View<bool*, HostMemorySpace> owned_host_;
};

} // namespace pcms
#endif // PCMS_MFEM_COUPLING_MFEM_ADAPTER_LAYOUT_H
