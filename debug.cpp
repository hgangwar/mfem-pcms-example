#include <Omega_h_mesh.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <gmsh.h>
#include <mpi.h>
#include <numeric>
#include <sstream>
#include <iostream>
#include <pcms/pcms.h>
/**
 * Save mesh elements to the separete meshes using rank2elem map
 * @param mesh        The full Omega_h mesh.
 * @param rank2elem   Mapping: rank → list of element IDs belonging to that rank.
 */
void save_ptn_mesh(Omega_h::Mesh& mesh,
                   const std::map<int, std::vector<int>>& rank2elem)
{
  const int nelems = mesh.nelems();
  Omega_h::Write<Omega_h::I32> elem_rank(nelems, -1);

  // Assign rank IDs to each element based on rank2elem mapping
  for (const auto& [rank, elems] : rank2elem) {
    for (auto eid : elems) {
      if (eid >= 0 && eid < nelems)
        elem_rank[eid] = rank;
      else
        fprintf(stderr, "[save_ptn_mesh] Warning: invalid element id %d for rank %d\n", eid, rank);
    }
  }

  // Remove existing "rank" tag if it exists
  if (mesh.has_tag(mesh.dim(), "rank"))
    mesh.remove_tag(mesh.dim(), "rank");

  // Add the new "rank" tag to elements (region dimension)
  mesh.add_tag<Omega_h::I32>(mesh.dim(), "rank", 1);
  mesh.set_tag<Omega_h::I32>(mesh.dim(), "rank", Omega_h::Read<Omega_h::I32>(elem_rank));

  printf("[save_ptn_mesh] Added rank tag to mesh elements.\n");

  // Save mesh in .osh format
  Omega_h::binary::write("partioned_cylinder.osh", &mesh);
}
void tagMeshElms(Omega_h::Mesh& mesh,
                      const pcms::Partition& partitions)
{
  const auto dim = mesh.dim();
  std::map<int, std::vector<int>> elemsPerRank;
  auto coords = mesh.coords();
  const auto elem2verts = mesh.ask_elem_verts();
  std::cout << "Number of elements: " << mesh.nelems() << std::endl;
  int counter = 0;
  for (int i = 0; i < mesh.nelems(); i++) {
    const auto tet_verts = Omega_h::gather_verts<4>(elem2verts, i);
    const auto tet_coords = Omega_h::gather_vectors<4,3>(coords, tet_verts);

    Omega_h::Vector<3> centroid = Omega_h::zero_vector<3>();
    for (int i = 0; i < 4; ++i) 
      centroid += tet_coords[i];
    centroid /= 4.0;
    std::array<double,3> points = {centroid[0], centroid[1], centroid[2]};
    
    auto rank = partitions.GetDr(i, dim, points);
    elemsPerRank[rank].push_back(i);
  }
  for ( const auto& [rank, elem_vec] : elemsPerRank ){
    printf("Element on rank %d : %d\n", rank, elem_vec.size());
    counter+=(elem_vec.size());
  }
  // Save the tagged mesh
  save_ptn_mesh(mesh, elemsPerRank);
  printf("Elements on all Ranks:%d\n", counter);
}
int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  MPI_Comm dup_comm;
  MPI_Comm_dup(MPI_COMM_WORLD, &dup_comm);

  Omega_h::Library lib(&argc, &argv, dup_comm);
  auto world = lib.world();
  int world_rank = world->rank();

  // ======= Binary::read =======
  {
    Omega_h::Mesh mesh(&lib);
    Omega_h::binary::read("../mesh/parmesh/cylinder.osh", world, &mesh);
    MPI_Barrier(world->get_impl());
    auto owned = mesh.owned(0);
    auto nents_owned = std::accumulate(owned.begin(), owned.end(), 0);
    std::stringstream ss;
    ss << "[binary::read] rank " << world_rank
       << " owned vertices = " << nents_owned << "\n";
    std::cout << ss.str();

    //Define Redev Partition
    redev::LO dim = 3; // 3D case
    redev::LOs ranks(2);
    std::iota(ranks.begin(),ranks.end(),0);
    redev::Reals cuts={0, 0.5};
    auto partition = redev::Partition{redev::RCBPtn{dim, ranks, cuts}};
    
    // do mesh migration to match the partition
    tagMeshElms(mesh, partition);
  }

  

  return 0;
}