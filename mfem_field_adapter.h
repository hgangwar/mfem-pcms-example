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
                    mfem::ParGridFunction& gf_data)
      : name_(std::move(name)),
        pmesh_(pmesh),
        pfes_(pfes),
        gf_data_(gf_data)
    {
      PCMS_FUNCTION_TIMER;
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
        // create a vector to store the serialized data
        mfem::Vector serialized_data(pfes_.GetTrueVSize());
        R->Mult(gf_data_, serialized_data);
      
        // ! instead of returning the serialized data, we need to write it to the buffer
        if (permutation.size() >0){ // check if permutation is empty
          for(int i=0; i<serialized_data.Size(); ++i){
            buffer[i] = serialized_data[permutation[i]];
          }
        }
        else {
          for(int i=0; i<serialized_data.Size(); ++i) {
            buffer[i] = serialized_data[i];
          }
        }
      }
      return pfes_.GetTrueVSize();
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

      // multiply the gf_data with the P matrix to get the deserialized data
      // create a vector to store the deserialized data
      //mfem::Vector deserialized_data(pfes_.GetTrueVSize());
      P->Mult(buffer_vector, gf_data_); // directly write to the gf_data_?
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
    auto * R = pfes_.GetRestrictionMatrix();
    if(!R) {
      std::cerr<<"R matrix is nullptr\n";
      std::abort();
    }
    mfem::Array<HYPRE_BigInt> gids;
    pmesh_.GetGlobalVertexIndices(gids);
    int size = gids.Size();
    mfem::Vector gid_vector(size);
    for(int i=0; i<size; ++i) {
      gid_vector[i] = gids[i];
    }
    mfem::Vector tgids(pfes_.GetTrueVSize());
    //R->BooleanMult(gids, tgids);
    R->Mult(gid_vector, tgids);
    //auto gids_host = gids.HostRead();
    // ? where should the return go?
    return {tgids.begin(), tgids.end()};
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
    // we need to create a counter for local index
    for (auto i = 0; i < vcoords.Size(); i+=3) { // class ids will be replaced with the node points
      std::array<double,3> coord{vcoords[i], vcoords[i+1], vcoords[i+2]};
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

private:
  std::string name_;
  mfem::ParMesh& pmesh_;
  mfem::ParFiniteElementSpace& pfes_;
  mfem::ParGridFunction& gf_data_;

};


} // namespace pcms

#endif // PCMS_COUPLING_XGC_FIELD_ADAPTER_H
