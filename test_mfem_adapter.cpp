#include "mfem_field_adapter.h"
#include <algorithm>
#include "mfem.hpp"
#include <vector>
#include <map>

using pcms::make_array_view;
using pcms::make_const_array_view;

/**
 * @brief Check if two values are close
 * @tparam T The type of the first value
 * @tparam T2 The type of the second value
 * @param val1 The first value
 * @param val2 The second value
 * @return True if the values are close, false otherwise
 * @note This is a dummy function
*/
template <typename T, typename T2>
bool is_close(T val1, T2 val2)
{
  // check if T2 is an integral type
  //static_assert(std::is_integral_v<T2>, "T2 should be counter/integral");
  // check if T is an integral type
  if constexpr (std::is_integral_v<T>) {
    return val1 == val2; // if T is an integral type, we can compare them with ==
  }
  return fabs(val1 - val2) < 1E-16; // if T is not an integral type, compare like floatig point comparisons
}

/**
 * @brief check if buffer and gf_data are close
 * @param buffer The buffer to check
 * @param gf_data The grid function data to check
 * @return True if the values are close, false otherwise
*/
bool check_data(const std::vector<double> & buffer, const mfem::Vector & vector_gf_data)
{
  // check if the size of the buffer is the same as the size of the gf_data
  if (buffer.size() != vector_gf_data.Size()) {
    return false;
  }
  // loop through the buffer and gf_data and check if the values are close
  for (int i = 0; i < buffer.size(); ++i) {
    if (!is_close(buffer[i], vector_gf_data[i])) {
      return false;
    }
  }
  return true;
}

/**
 * @brief Create a true grid function data
 * @param gf_data The grid function data
 * @return The true grid function data
*/
mfem::Vector make_true_gf_data(const mfem::GridFunction & gf_data, const mfem::ParFiniteElementSpace & pfes)
{
  // create a vector to store the true grid function data
  mfem::Vector true_gf_data(gf_data.Size());
  // get restriction matrix
  auto * R = pfes.GetRestrictionMatrix();
  R->Mult(gf_data, true_gf_data);
  return true_gf_data;
}


int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  // construct/load a parallel mfem mesh
	mfem::Mesh mesh("../mesh/testmesh.msh");
	//mfem::ParMesh *pmesh = new mfem::ParMesh(MPI_COMM_WORLD, *mesh);
  mfem::ParMesh pmesh(MPI_COMM_WORLD, mesh);
  //delete mesh;
  mesh.Clear();
	//
	//
	//
  // load a finite element space
  // for fes: we need the mesh and the finite element collection
  const int dim = 3; // 3D case
  int order = 1;
  //mfem::FiniteElementCollection *fec;
  //fec = new mfem::H1_FECollection(order, dim);
  mfem::H1_FECollection fec(order, dim);

	//mfem::ParFiniteElementSpace *pfes = new mfem::ParFiniteElementSpace(pmesh, fec);
  mfem::ParFiniteElementSpace pfes(&pmesh, &fec);
	
  // auto pfes = ...
  // construct some dummy data on pfes finite element space


	mfem::ParGridFunction gf_data(&pfes);

	// we have to create some coefficients to project on the gf_data
	//mfem::ConstantCoefficient three(3.0); 
	//gf_data.ProjectCoefficient(three);
	// TODO: Will also do function coefficient
	std::iota(gf_data.begin(), gf_data.end(), 0.0);
	
  // gf_data is a "vector" fill it up with something we know
  // std::iota(gf_data.begin(), gf_data.end(), 0.0); // 
  // ! we don't need this since we can fill gf_data with projection function
  // construct mfem field adapter
  pcms::MFEMFieldAdapter adapter(std::string("mfem_field_adapter"), pmesh, pfes, gf_data);
  // want to call serialize
  std::vector<double> buffer;
  std::vector<int> permutation;
  auto size = adapter.Serialize(make_array_view(buffer), make_const_array_view(permutation));
  // this isn't quite right yet...need the TDOF size
  auto tsize = pfes.GetTrueVSize();
  PCMS_ALWAYS_ASSERT(size == tsize);
  buffer.resize(size); 

  // should copy data from gf_data into buffer vector
  // note we pass empty permutatation array, so adapter needs to handle this case
  adapter.Serialize(make_array_view(buffer), make_const_array_view(permutation));
  
  // check data function returns true if all entries in buffer are close to
  // those in the gf_data
  auto true_gf_data = make_true_gf_data(gf_data, pfes);
  PCMS_ALWAYS_ASSERT(check_data(buffer, true_gf_data));

  // modify values in our buffer
  std::transform(buffer.begin(), buffer.end(), buffer.begin(), [](double val){ return val*2; });
  adapter.Deserialize(make_const_array_view(buffer), make_const_array_view(permutation));
  true_gf_data = make_true_gf_data(gf_data, pfes);
  PCMS_ALWAYS_ASSERT(check_data(buffer, gf_data));

  // check that GIDS is correct
  // get gids directly from the mesh and compare values to GetGids function
  // get gids from mesh, loop through gids from mesh and from function and make sure they are equal
  auto *R = pfes.GetRestrictionMatrix();
  mfem::Array<HYPRE_BigInt>gids;
  pmesh.GetGlobalVertexIndices(gids);
  int size_gids = gids.Size();
  mfem::Vector gid_vector(size_gids);
  for (int i = 0; i < size_gids ; ++i) {
     gid_vector[i] = gids[i];
  } 
  mfem::Vector tgids(pfes.GetTrueVSize());
  R->Mult(gid_vector,tgids);
  auto gids_adapter = adapter.GetGids();
  //int gIndx; 
  PCMS_ALWAYS_ASSERT(tgids.Size() == gids_adapter.size());
  for (int i = 0 ; i < gids_adapter.size(); ++i) {
   PCMS_ALWAYS_ASSERT(gids_adapter[i] == tgids[i]);
     
  }
  
  // * this section will create a partition and check if the reverse partition map is correct
  // create an RCB partition
  
  std::vector<pcms::LO> ranks(8);
  std::iota(ranks.begin(),ranks.end(),0);
  std::vector<pcms::Real> cuts = {0,/*x*/0.5,/*y*/0.75,0.25,/*z*/0.1,0.4,0.8,0.3};
  auto partition = redev::Partition{redev::RCBPtn{dim, ranks, cuts}}; // ? partition constructor work here? How?
  auto& rcb_partition = std::get<redev::RCBPtn>(partition);
  // check that reverse partition map holds the data we expect
  // create a partition in our case should be RCB
  auto reverse_partition_map = adapter.GetReversePartitionMap(partition);
  
  // ? Now? What to do? how do I know that which vertex is in which rank?
  // based on the values for the partition that we constructed verify that the
  // entries in reverse partition map are sending to the correct rank 
   int local_index = 0;
   mfem::Vector vcoords;
   pmesh.GetVertices(vcoords);   
   for (auto i = 0; i < vcoords.Size(); i+=3){
       auto coord = std::array<double, 3>{vcoords[i], vcoords[i+1], vcoords[i+2]};
   	int to_rank = rcb_partition.GetRank(coord);
   	// check that the reverse parition map has key that is to_rank
   	// this if is checking if the to_rank rank is present in the list of keys of
   	// reverse partition map
   	//
   	if(auto it = reverse_partition_map.find(to_rank); it != reverse_partition_map.end()) {
          const auto& local_indices = it->second;
          auto local_index_it = std::find(local_indices.begin(), local_indices.end(), local_index);
          if(local_index_it == local_indices.end()) {
            std::cerr<<"Reverse partition map for sending to rank "<<to_rank<<" is missing local index "<<local_index<<"\n";
            std::cerr<<"reverse partition map contains:["; 
          for(auto v : local_indices) {
            std::cerr<<v<<" ";
          }
          std::cerr<<"]\n";
          return 1;
        }
   	}
   	else {
   	  std::cerr<<"to rank key " << to_rank <<" does not exist\n";
          std::cerr<<"Current keys are: ";
          for(const auto& vls: reverse_partition_map) {
            std::cerr<<vls.first << " ";
          }
          std::cerr<<"\n";
   	  return 1;
   	}
   	++local_index;
  }
  
MPI_Finalize();
}
