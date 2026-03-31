#ifndef PCMS_COUPLING_MFEM_FIELD_ADAPTER_H
#define PCMS_COUPLING_MFEM_FIELD_ADAPTER_H
#include <pcms/utility/types.h>
#include <pcms/utility/memory_spaces.h>
#include <pcms/field.h>
#include <vector>
#include <redev_variant_tools.h>
#include <pcms/utility/assert.h>
#include <pcms/utility/array_mask.h>
#include <pcms/utility/profile.h>
#include <pcms/partition.h>
#include "mfem.hpp"

namespace pcms
{
  class MFEMFieldAdapter
  {
  public:
    using memory_space = HostMemorySpace;
    using value_type = Real;
    using coordinate_element_type = Real;
    /**
    *
    * @param name name of the field
    * @param gf_data the mfem grid function data
    */
    MFEMFieldAdapter(std::string name,
                 mfem::ParMesh& pmesh,
                 mfem::ParFiniteElementSpace& pfes,
                 mfem::ParGridFunction& gf_data,
                 bool use_mask = false,
                 pcms::LO attr = -1)
  : name_(std::move(name)),
    pmesh_(pmesh),
    pfes_(pfes),
    gf_data_(gf_data)
    {
      PCMS_FUNCTION_TIMER;

      if (use_mask)
      {
        PCMS_ALWAYS_ASSERT(attr >= 0);
        create_mask(attr);
      }
      else
      {
        hasmask_ = false;
        packed_size_ = 0;

        mask_storage_.SetSize(0);
        mask_view_ = Rank1View<pcms::LO, Kokkos::HostSpace>(nullptr, 0);
      }
    }
    bool create_mask(pcms::LO attr)
    {
      PCMS_FUNCTION_TIMER;

      const int ne = pmesh_.GetNE();
      const int nv = pmesh_.GetNV();

      printf("Number of elements in the mesh: %d\n", ne);
      printf("Number of vertices in the mesh: %d\n", nv);

      mask_storage_.SetSize(nv);
      mask_storage_ = 0;

      pcms::LO count = 0;
      mfem::Array<int> vert_ids;

      for (int e = 0; e < ne; ++e)
      {
        if (pmesh_.GetAttribute(e) != attr) { continue; }

        pmesh_.GetElementVertices(e, vert_ids);

        for (int j = 0; j < vert_ids.Size(); ++j)
        {
          const int v = vert_ids[j];

          if (mask_storage_[v] == 0)
          {
            mask_storage_[v] = ++count;
          }
        }
      }

      packed_size_ = count;
      hasmask_ = (packed_size_ > 0);

      printf("Filtered %d unique vertices.\n", packed_size_);

      mask_view_ = Rank1View<pcms::LO, Kokkos::HostSpace>(
          mask_storage_.GetData(), mask_storage_.Size());

      return hasmask_;
    }
   
    // REQUIRED
    // serialize the data
    int Serialize(
      Rank1View<value_type, memory_space> buffer,
      Rank1View<const pcms::LO, memory_space> permutation) const
    {
      PCMS_FUNCTION_TIMER;
      static_assert(std::is_same_v<memory_space, pcms::HostMemorySpace>,
                    "gpu space unhandled\n");
      if(buffer.size() >0) {
        // get restriction matrix
        auto * R = pfes_.GetRestrictionMatrix();
        if(!R) {
          std::cerr<<"R matrix is nullptr\n";
          std::abort();
        } 
        // multiply the gf_data with the R matrix to get the serialized data
        mfem::Vector serialized_data(pfes_.GetTrueVSize());
        R->Mult(gf_data_, serialized_data);
        pcms::LO filtered_size = has_mask()?packed_size_:pfes_.GetTrueVSize();
        mfem::Vector filtered_data(filtered_size);
        printf("\n Size of filtered_data : %d, mask size: %d\n", filtered_data.Size(), this->mask_storage_.Size());

        if (has_mask()) {
          for (pcms::LO i = 0; i < serialized_data.Size(); ++i) {
            if (mask_view_(i) > 0) {
              const pcms::LO idx = mask_view_(i) - 1;

              if (idx >= packed_size_) {
                std::cerr << "Assertion failed: index out of range for filtered array. "
                          <<" itr:"<<i<< " idx: " << idx << ", packed_size_: " << packed_size_ << std::endl;
                std::abort();
              }
              filtered_data[idx] = serialized_data[i];
            }
          }
        } else {
          for (pcms::LO i = 0; i < serialized_data.Size(); ++i) {
            filtered_data[i] = serialized_data[i];
          }
        }
        // ! instead of returning the serialized data, we need to write it to the buffer
        if (permutation.size() >0){ // check if permutation is empty
          for(int i=0; i<filtered_data.Size(); ++i){
            buffer[i] = filtered_data[permutation[i]];
          }
        }
        else {
          for(int i=0; i<filtered_data.Size(); ++i) {
            buffer[i] = filtered_data[i];
          }
        }
      }

      return packed_size_;
    }

    // REQUIRED
    // deserialize the data
    void Deserialize(
      Rank1View<const value_type, memory_space> buffer,
      Rank1View<const pcms::LO, memory_space> permutation) const
    {
      PCMS_FUNCTION_TIMER;
      static_assert(std::is_same_v<memory_space, pcms::HostMemorySpace>,
                    "gpu space unhandled\n");
      //if (RankParticipatesCouplingCommunication()) {
      //  // ? just need to replace data_ with the powerdensity_
      //  mask_.ToFullArray(buffer, gf_data_, permutation);
      //}

      // get the prolongation matrix
      auto * P = pfes_.GetProlongationMatrix();
      if(!P) {
              std::cerr<<"P matrix is nullptr\n";
              std::abort();
      }
      mfem::Vector buffer_vector(buffer.size());
      if (permutation.size() >0) {
        for(int i=0; i<buffer.size(); ++i) {
          buffer_vector[i] = buffer[permutation[i]];
        }
      }
      else {
        for(int i=0; i<buffer.size(); ++i) {
          buffer_vector[i] = buffer[i];
        }
      }
      auto * R = pfes_.GetRestrictionMatrix();
      mfem::Vector serialized_data(pfes_.GetTrueVSize());
      R->Mult(gf_data_, serialized_data);
      if (has_mask()) {
        // Merge the buffer and true data

      }
      for (int i=0; i<serialized_data.Size(); ++i) {
        if (!has_mask()) {
          serialized_data[i] = buffer_vector[i];
        }
        else if (mask_view_(i)) {
          const pcms::LO idx = mask_view_(i) - 1;
          serialized_data[i] = buffer_vector[idx];
        }
      }

      // multiply the gf_data with the P matrix to get the deserialized data
      // create a vector to store the deserialized data
      //mfem::Vector deserialized_data(pfes_.GetTrueVSize());
      P->Mult(serialized_data, gf_data_); // directly write to the gf_data_?
      // Synchronize fields? Think mult with P will handle parallel synchronization of the owned/ghost DOF data
    }

 /**
  * @brief Get the Gids
  * @return std::vector<GO>
  * 
 */
    [[nodiscard]] std::vector<GO> GetGids() const
    {
      PCMS_FUNCTION_TIMER;
      mfem::Array<HYPRE_BigInt> gids;
      pmesh_.GetGlobalVertexIndices(gids);

      if (has_mask()) {
        std::vector<GO> filtered_gids(packed_size_);
        for (int i = 0; i < mask_storage_.Size(); ++i) {
          if (mask_view_(i) > 0) {
            filtered_gids[mask_view_(i) - 1] = static_cast<GO>(gids[i]);
          }
        }
        return filtered_gids;
      }
      return {gids.begin(), gids.end()};
    }
  // REQUIRED
  [[nodiscard]] ReversePartitionMap GetReversePartitionMap(
    const pcms::Partition& partition) const
  {
    PCMS_FUNCTION_TIMER;
    pcms::ReversePartitionMap reverse_partition;
   // note GetVertices assumes that the mesh is not higher order
   // if we have a higher order mesh, we need to use GetNodes 
    mfem::Vector vcoords;
    pcms::LO dim = pmesh_.Dimension();
    pmesh_.GetVertices(vcoords);
    int local_index=0;
    std::array<double, 3> coord;
    // we need to create a counter for local index
    for (auto i = 0; i < vcoords.Size(); i+=dim) { // class ids will be replaced with the node points
      //std::array<double,3> coord{vcoords[i], vcoords[i+1], vcoords[i+2]};
      std::copy(vcoords.begin() + i, vcoords.begin() + i + dim, coord.begin());

      auto dr = partition.GetDr(local_index, dim, coord);
      reverse_partition[dr].emplace_back(local_index++); // it should be some counter since it is going 3 at a time
    }     
    int counter = 0;
    for ( const auto& [rank, vertx] : reverse_partition ){
      printf("Vertex in Reverse Partition map, rank %d : %d\n", rank, vertx.size());
      counter+=(vertx.size());
    }

  return reverse_partition; 
  }

    const Rank1View<pcms::LO, Kokkos::HostSpace>& get_mask_view() const
    {
      return mask_view_;
    }

    pcms::LO get_packed_size() const
    {
      return packed_size_;
    }

    bool has_mask() const
    {
      return hasmask_;
    }
private:
  std::string name_;
  mfem::ParMesh& pmesh_;
  mfem::ParFiniteElementSpace& pfes_;
  mfem::ParGridFunction& gf_data_;
  mfem::Array<pcms::LO> mask_storage_;
  Rank1View<pcms::LO, Kokkos::HostSpace> mask_view_;

  bool hasmask_ = false;
  pcms::LO packed_size_ = 0;

};


} // namespace pcms

#endif // PCMS_COUPLING_XGC_FIELD_ADAPTER_H
