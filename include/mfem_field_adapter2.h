#ifndef PCMS_MFEM_COUPLING_MFEM_FIELD_ADAPTER2_H
#define PCMS_MFEM_COUPLING_MFEM_FIELD_ADAPTER2_H

#include <pcms/field.h>
#include <pcms/utility/assert.h>
#include <pcms/utility/memory_spaces.h>
#include <pcms/utility/profile.h>
#include <pcms/utility/types.h>

#include "mfem_adapter_layout.h"

#include <mfem.hpp>

namespace pcms
{

template <typename T>
class MFEMFieldAdapter2 : public FieldT<T>
{
public:
  using memory_space = HostMemorySpace;
  using value_type = T;

  explicit MFEMFieldAdapter2(const MFEMFieldsAdapterLayout& layout);

  const FieldLayout& GetLayout() const override;

  int Serialize(
    Rank1View<T, HostMemorySpace> buffer,
    Rank1View<const LO, HostMemorySpace> permutation) const override;

  void Deserialize(Rank1View<const T, HostMemorySpace> buffer,
                   Rank1View<const LO, HostMemorySpace> permutation) override;

  Rank1View<const T, HostMemorySpace> GetDOFHolderData() const override;

  void SetDOFHolderData(Rank1View<const T, HostMemorySpace> data) override;

  LocalizationHint GetLocalizationHint(
    CoordinateView<HostMemorySpace>) const override
  {
    throw pcms_error("MFEMFieldAdapter2 is communication-only.");
  }

  void Evaluate(LocalizationHint,
                FieldDataView<T, HostMemorySpace>) const override
  {
    throw pcms_error("MFEMFieldAdapter2 is communication-only.");
  }

  void EvaluateGradient(FieldDataView<T, HostMemorySpace>) override
  {
    throw pcms_error("MFEMFieldAdapter2 is communication-only.");
  }

  bool CanEvaluateGradient() override { return false; }

private:
  void AssertVertexOnlyField() const;

private:
  const MFEMFieldsAdapterLayout& layout_;

  mfem::ParMesh& pmesh_;
  mfem::ParFiniteElementSpace& pfes_;
  mfem::ParGridFunction& gf_;

  mutable Kokkos::View<T*, HostMemorySpace> data_cache_;
};

template <typename T>
MFEMFieldAdapter2<T>::MFEMFieldAdapter2(const MFEMFieldsAdapterLayout& layout)
  : layout_(layout),
    pmesh_(layout.GetMesh()),
    pfes_(layout.GetFESpace()),
    gf_(layout.GetGridFunction())
{
  PCMS_FUNCTION_TIMER;

  AssertVertexOnlyField();

  data_cache_ = Kokkos::View<T*, HostMemorySpace>(
    "mfem_field_adapter_cache", layout_.GetNumOwnedDofHolder());
}

template <typename T>
void MFEMFieldAdapter2<T>::AssertVertexOnlyField() const
{
  const auto nodes = layout_.GetNodesPerDim();

  PCMS_ALWAYS_ASSERT(nodes[0] == 1);
  PCMS_ALWAYS_ASSERT(nodes[1] == 0);
  PCMS_ALWAYS_ASSERT(nodes[2] == 0);
  PCMS_ALWAYS_ASSERT(nodes[3] == 0);

  PCMS_ALWAYS_ASSERT(layout_.GetNumComponents() == 1);
  PCMS_ALWAYS_ASSERT(pfes_.GetVDim() == 1);
}

template <typename T>
const FieldLayout& MFEMFieldAdapter2<T>::GetLayout() const
{
  return layout_;
}

template <typename T>
Rank1View<const T, HostMemorySpace> MFEMFieldAdapter2<T>::GetDOFHolderData()
  const
{
  PCMS_FUNCTION_TIMER;

  mfem::Array<int> vdofs;

  if (layout_.HasMask()) {
    const auto mask = layout_.GetMask();

    for (int v = 0; v < pmesh_.GetNV(); ++v) {
      if (mask[v] == 0) {
        continue;
      }

      pfes_.GetVertexDofs(v, vdofs);
      PCMS_ALWAYS_ASSERT(vdofs.Size() == 1);

      const LO out = mask[v] - 1;
      const int dof = vdofs[0];

      data_cache_(out) = static_cast<T>(gf_(dof));
    }
  } else {
    for (int v = 0; v < pmesh_.GetNV(); ++v) {
      pfes_.GetVertexDofs(v, vdofs);
      PCMS_ALWAYS_ASSERT(vdofs.Size() == 1);

      const int dof = vdofs[0];

      data_cache_(v) = static_cast<T>(gf_(dof));
    }
  }

  return make_const_array_view(data_cache_);
}

template <typename T>
void MFEMFieldAdapter2<T>::SetDOFHolderData(
  Rank1View<const T, HostMemorySpace> data)
{
  PCMS_FUNCTION_TIMER;

  PCMS_ALWAYS_ASSERT(static_cast<LO>(data.size()) ==
                     layout_.GetNumOwnedDofHolder());

  mfem::Array<int> vdofs;

  if (layout_.HasMask()) {
    const auto mask = layout_.GetMask();

    for (int v = 0; v < pmesh_.GetNV(); ++v) {
      if (mask[v] == 0) {
        continue;
      }

      pfes_.GetVertexDofs(v, vdofs);
      PCMS_ALWAYS_ASSERT(vdofs.Size() == 1);

      const LO in = mask[v] - 1;
      const int dof = vdofs[0];

      gf_(dof) = static_cast<double>(data[in]);
    }
  } else {
    for (int v = 0; v < pmesh_.GetNV(); ++v) {
      pfes_.GetVertexDofs(v, vdofs);
      PCMS_ALWAYS_ASSERT(vdofs.Size() == 1);

      const int dof = vdofs[0];

      gf_(dof) = static_cast<double>(data[v]);
    }
  }
}

template <typename T>
int MFEMFieldAdapter2<T>::Serialize(
  Rank1View<T, HostMemorySpace> buffer,
  Rank1View<const LO, HostMemorySpace> permutation) const
{
  PCMS_FUNCTION_TIMER;

  const auto data = GetDOFHolderData();

  PCMS_ALWAYS_ASSERT(buffer.size() == data.size());

  for (size_t i = 0; i < data.size(); ++i) {
    const LO dst = permutation.size() > 0 ? permutation[i] : i;
    buffer[dst] = data[i];
  }

  return static_cast<int>(data.size());
}

template <typename T>
void MFEMFieldAdapter2<T>::Deserialize(
  Rank1View<const T, HostMemorySpace> buffer,
  Rank1View<const LO, HostMemorySpace> permutation)
{
  PCMS_FUNCTION_TIMER;

  PCMS_ALWAYS_ASSERT(static_cast<LO>(buffer.size()) ==
                     layout_.GetNumOwnedDofHolder());

  if (permutation.size() == 0) {
    SetDOFHolderData(buffer);
    return;
  }

  Kokkos::View<T*, HostMemorySpace> sorted("mfem_sorted_buffer", buffer.size());

  for (size_t i = 0; i < buffer.size(); ++i) {
    sorted(i) = buffer[permutation[i]];
  }

  SetDOFHolderData(make_const_array_view(sorted));
}

} // namespace pcms

#endif