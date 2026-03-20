#include <iostream>
#include <Omega_h_mesh.hpp>
#include <Omega_h_file.hpp>
#include <gmsh.h>
#include <numeric>
#include <sstream>


int main(int argc, char** argv) {

  Omega_h::Library lib(&argc, &argv);
  auto world = lib.world();              
  MPI_Comm mpi_comm = world->get_impl();
  int world_rank;
  MPI_Comm_rank(mpi_comm, &world_rank);

  // ======= Case 1: binary::read =======
  {
    Omega_h::Mesh mesh(&lib);
    Omega_h::binary::read("cylinder.osh", world, &mesh);

    MPI_Barrier(mpi_comm);
    auto owned = mesh.owned(0);
    auto nents_owned = std::accumulate(owned.begin(), owned.end(), 0);
    std::stringstream ss;
    ss << "[binary::read] rank " << world_rank
       << " owned vertices = " << nents_owned << "\n";
    std::cout << ss.str();
  }

  // ======= Case 2: gmsh::read_parallel =======
  {
    gmsh::initialize(argc, argv, false);
    auto mesh = Omega_h::gmsh::read("cylinder.msh", world);

    MPI_Barrier(mpi_comm);
    auto owned = mesh.owned(0);
    auto nents_owned = std::accumulate(owned.begin(), owned.end(), 0);
    std::stringstream ss;
    ss << "[gmsh::read] rank " << world_rank
       << " owned vertices = " << nents_owned << "\n";
    std::cout << ss.str();

    gmsh::finalize();
  }
  return 0;
}

