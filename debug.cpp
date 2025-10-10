#include <Omega_h_mesh.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <gmsh.h>
#include <mpi.h>
#include <numeric>
#include <sstream>
#include <iostream>

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  MPI_Comm dup_comm;
  MPI_Comm_dup(MPI_COMM_WORLD, &dup_comm);

  Omega_h::Library lib(&argc, &argv, dup_comm);
  auto world = lib.world();
  int world_rank = world->rank();

  // ======= Case 1: binary::read =======
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
  }

  
  // ======= Case 2: gmsh::read_parallel =======
  {
    gmsh::initialize();
    auto mesh = Omega_h::gmsh::read("../mesh/parmesh/cylinder.msh", world);
    MPI_Barrier(world->get_impl());
    auto owned = mesh.owned(0); 
    auto nents_owned = std::accumulate(owned.begin(), owned.end(), 0);
    std::stringstream ss;
    ss << "[gmsh::read] rank " << world_rank
       << " owned vertices = " << nents_owned << "\n";
    std::cout << ss.str();
    gmsh::finalize();
  }

  // Clean up
  MPI_Comm_free(&dup_comm);
  MPI_Finalize();

  return 0;
}