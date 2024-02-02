/* ! Description:
 * This file is a test for the coupling between two MFEM solvers using PCMS.
 * Both of the solvers are based on the same mesh.
 * The thermal solver just solves simple heat equation. And using the solution
 * of the thermal solver (density using a eof) the flux solver solves the
 * diffusion equation of neutron flux.
 * The coupling is done by PCMS.
*/

/* ! Test description
 * This test is based on Catch2 framework.
 * The test is done by creating two dummy solvers (thermal and flux) and
 * coupling them using PCMS. The sent and received data are compared to
 * check the correctness of the coupling.
 * The test is now done for 3D meshes.
*/
#include <iostream>
#include <cmath>
#include <mpi.h>
//#include <catch2/catch_session.hpp>
//#include <catch2/catch_test_macros.hpp>
#include "mfem.hpp"
#include "mfem_field_adapter.h"
#include <Omega_h_mesh.hpp>
#include <pcms/pcms.h>
#include <pcms/types.h>
#include <Omega_h_file.hpp>
#include <Omega_h_for.hpp>
#include <redev_variant_tools.h>
#include <pcms/omega_h_field.h>
#include <gmsh.h>
#include <gmsh.h>
#include "Omega_h_element.hpp"
#include "Omega_h_shape.hpp"
#include <sstream>
#include "test_support.h"
//#include <atomic>

//using namespace Omega_h;

using pcms::Copy;
using pcms::CouplerClient;
using pcms::CouplerServer;
using pcms::FieldEvaluationMethod;
using pcms::FieldTransferMethod;
using pcms::GO;
using pcms::Lagrange;
using pcms::make_array_view;
using pcms::OmegaHFieldAdapter;
using pcms::MFEMFieldAdapter;

// this function checks if two doubles are close to each other
// within a tolerance

bool is_close(double a, double b, double tol = 1e-12)
{
  return std::abs(a - b) < tol;
}

/* @brief function to get number of ranks in a communicator object
  * @param[in] comm: MPI communicator
  * @return number of ranks in the communicator
*/
int get_comm_size(MPI_Comm comm)
{
  int comm_size;
  MPI_Comm_size(comm, &comm_size);
  return comm_size;
}

/* @brief function to create cut positions for number of ranks in 3D
  * @param[in] num_ranks: number of ranks in the communicator
  * @return vector of cut positions
  * @details the cut positions works like this:
  *         for only one rank, the cut position vector is {0}
  *         for two ranks, the cut position vector is {0, 0.5}
  *         for four ranks, the cut position vector is {0, 0.5, 0.75, 0.25}
  *         for eight ranks, the cut position vector is {0, 0.5, 0.75, 0.25, 0.1, 0.4, 0.8, 0.3}
  *         and so on
  *        the cut positions are in the range of [0,1]
*/
std::vector<pcms::Real> get_cut_positions(int num_ranks)
{
  // this function only works for 3D meshes
  // for now, it only works upto 8 ranks, so check the number of ranks
  if (num_ranks > 8)
  {
    std::cout << "The number of ranks is more than 8\n";
    std::cout << "The get_cut_position function only works for 8 ranks or less\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  // create a vector of cut positions of size num_ranks
  std::vector<pcms::Real> cuts(num_ranks);
  // fill the first cut position with 0
  cuts[0] = 0;
  
  // fill the rest of the cut positions with 0.5
  if (num_ranks > 1){
    for (int i = 1; i < num_ranks; i++)
    {
      cuts[i] = 0.5;
    }
  }
  return cuts;
}


void mfem_coupler(MPI_Comm comm, Omega_h::Mesh& mesh)
{
  // create partition to give it to the server coupler class constructor
  const int dim = 3; // 3D case
  //std::vector<pcms::LO> ranks(8); // does it need to match with # of ranks of comm

  int comm_size = get_comm_size(comm);
  std::cout << "The number of ranks is: " << comm_size << "\n";
  std::vector<pcms::LO> ranks(comm_size);
  std::iota(ranks.begin(),ranks.end(),0);

  // check if the comm_size is power of 2
  if ((comm_size & (comm_size - 1)) != 0)
  {
    std::cout << "The number of ranks supplied to the coupler is not power of 2\n";
    std::cout << "In needs power of 2 since Omega_h follows tree partition scheme\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  //std::vector<pcms::Real> cuts = {0,/*x*/0.5,/*y*/0.75,0.25,/*z*/0.1,0.4,0.8,0.3};
  //std::vector<pcms::Real> cuts = {0};
  std::vector<pcms::Real> cuts = get_cut_positions(comm_size);

  auto partition = redev::Partition{redev::RCBPtn{dim, ranks, cuts}};
  auto& rcb_partition = std::get<redev::RCBPtn>(partition);

  // create the recursive partition
  test_support::RecursivePartition recursive_partition;
  recursive_partition.ranks = ranks;
  recursive_partition.cuts = cuts;
  // do mesh migration to match the partition
  test_support::migrateMeshElms(mesh, recursive_partition);

  pcms::CouplerServer cpl("mfem_couple_server", comm, partition, mesh); 

  auto* flux_app = cpl.AddApplication("fluxClient");
  auto* thermal_app = cpl.AddApplication("thermalClient");
  
  std::cout << "coupler: coupler server and applications are created \n";

  // is_overlap is a vector of size mesh.nents(0) and is initialized to 1
  Omega_h::Write<Omega_h::I8> is_overlap(mesh.nents(0));
  Omega_h::parallel_for(is_overlap.size(), OMEGA_H_LAMBDA(int i)
  {
    is_overlap[i] = 1;
  });


  auto* flux_density_field = flux_app->AddField(
    "density", OmegaHFieldAdapter<double>("flux_density", mesh, is_overlap),
    FieldTransferMethod::Copy, // to Omega_h
    FieldEvaluationMethod::None,
    FieldTransferMethod::Copy, // from Omega_h
    FieldEvaluationMethod::None, is_overlap);

  auto* thermal_density_field = thermal_app->AddField(
    "density", OmegaHFieldAdapter<double>("flux_density", mesh, is_overlap),
    FieldTransferMethod::Copy, FieldEvaluationMethod::None,
    FieldTransferMethod::Copy, FieldEvaluationMethod::None, is_overlap);
  
  std::cout << "coupler: fields are created \n";

  std::cout << "coupler: starting receive from thermalSolver\n";
  thermal_app->ReceivePhase([&]() { thermal_density_field->Receive(); });
  std::cout << "coupler: done receiving, starting send\n";
  flux_app->SendPhase([&]() { flux_density_field->Send(); });
  std::cout << "coupler: done sending to fluxSolver\n";
}

// thermal solver function for dummy thermal solver
// this function just loads the mesh in parallel and
// adds a dummy grid function to the mesh and sends the data
// to the flux solver using PCMS
//
// @param[in] mesh_file_name: the name of the mesh file
// @param[in] comm: MPI communicator
// TODO: instead to copying the data from gf to pgf
//      we can return the pointer to the grid function
//      which will save the time of copying the data
int thermal_solver(const std::string& mesh_file_name, MPI_Comm comm)
{
  // load the mesh using mfem
  mfem::Mesh *mesh = new mfem::Mesh(mesh_file_name.c_str(), 1, 1);
  mfem::ParMesh *pmesh = new mfem::ParMesh(comm, *mesh);
  int dim = pmesh->Dimension();

  // create the fe space
  mfem::FiniteElementCollection *fec = new mfem::H1_FECollection(1, dim);
  mfem::ParFiniteElementSpace *fespace = new mfem::ParFiniteElementSpace(pmesh, fec);

  // create the grid function
  mfem::GridFunction gf(fespace); 
  // fill it with 1.23 value using iota
  std::iota(gf.begin(), gf.end(), 1.23);
  // make it a parallel grid function
  mfem::ParGridFunction pgf = mfem::ParGridFunction(fespace, gf);


  // create the PCMS client
  CouplerClient cpl("thermalClient", comm);
  cpl.AddField("density", MFEMFieldAdapter(std::string("thermal_density"), *pmesh, *fespace, pgf));
  cpl.BeginSendPhase();
  cpl.SendField("density");
  cpl.EndSendPhase();



  // delete the mesh and fe space
  delete fespace;
  delete fec;
  delete pmesh;
  delete mesh;
  return 0;
}

// flux solver function for dummy flux solver
// this function is exactly the same as the thermal solver
// except that it receives the data from the thermal solver
// using PCMS
// @param[in] mesh_file_name: the name of the mesh file
// @param[in] comm: MPI communicator
//
// @return the grid function of the flux solver

int flux_solver(const std::string& mesh_file_name, MPI_Comm comm)
{
  // load the mesh using mfem
  mfem::Mesh *mesh = new mfem::Mesh(mesh_file_name.c_str(), 1, 1);
  mfem::ParMesh *pmesh = new mfem::ParMesh(comm, *mesh);
  int dim = pmesh->Dimension();

  // create the fe space
  mfem::FiniteElementCollection *fec = new mfem::H1_FECollection(1, dim);
  mfem::ParFiniteElementSpace *fespace = new mfem::ParFiniteElementSpace(pmesh, fec);

  // create the reference grid function
  mfem::GridFunction ref_gf(fespace); 
  // fill it with 1.23 value using iota
  std::iota(ref_gf.begin(), ref_gf.end(), 1.23);

  // this time we receive the data from the thermal solver
  // so the grid function is not filled
  // make it a parallel grid function
  mfem::ParGridFunction pgf = mfem::ParGridFunction(fespace);

  // create the PCMS client
  CouplerClient cpl("fluxClient", comm);
  cpl.AddField("density", MFEMFieldAdapter(std::string("flux_density"), *pmesh, *fespace, pgf));
  cpl.BeginReceivePhase();
  cpl.ReceiveField("density");
  cpl.EndReceivePhase();

  // copy the data from pgf to gf
  mfem::GridFunction gf(pgf);

  // compare the size of the grid functions
  if (gf.Size() != ref_gf.Size())
  {
    std::cout << "The size of the grid functions are not the same\n";
    std::cout << "The size of the sent grid function is: " << ref_gf.Size() << "\n";
    std::cout << "The size of the grid function is: " << gf.Size() << "\n";
    std::cout << "The test failed\n";
    return 1;
  }
  if (gf.Size() == ref_gf.Size())
  {
    std::cout << "The size of the grid functions are the same\n";
  }

  int flux_procs;
  MPI_Comm_size(comm, &flux_procs);
  int flux_myid;
  MPI_Comm_rank(comm, &flux_myid);

  if (flux_myid == 0)
  {
    for (int i = 0; i < gf.Size(); i++)
    {
      if (!is_close(gf[i], ref_gf[i]))
      {
        std::cout << "The data of the grid functions are not the same\n";
        std::cout << "The test failed\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
    }
  }
  // delete the mesh and fe space
  delete fespace;
  delete fec;
  delete pmesh;
  delete mesh;
  return 0;
}


// the coupler_run function loads the mesh using Omega_h and
// calls the thermal and flux solvers
// and the coupler server and catches the grid functions
// then compares the data

int main(int argc, char *argv[])
{
  // initialize MPI
  MPI_Init(&argc, &argv);
  int world_rank, world_size;

  // get the world rank and size
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // if the world size is not 3 then return 1
  //if (world_size != 3)
  //{
  //  std::cout << "World size must be 3 for " << argv[0] << std::endl;
  //  MPI_Abort(MPI_COMM_WORLD, 1);
  //}

  // number of ranks for each solver
  int thermal_solver_ranks = 1;
  int flux_solver_ranks = 1;
  int coupler_server_ranks = world_size - thermal_solver_ranks - flux_solver_ranks;
  // reduce the coupler server ranks to the nearest power of 2
  coupler_server_ranks = std::pow(2, std::floor(std::log2(coupler_server_ranks)));

  // split the communicator based on the rank of each process
  //MPI_Comm comm;
  //int color = world_rank / 3; // Determine color 
  //MPI_Comm_split(MPI_COMM_WORLD, color, world_rank, &comm);


  // mesh name
  std::string mesh_file_name = "../mesh/parmesh/cylinder.msh";
  
  MPI_Comm server_comm, flux_comm, thermal_comm;
  MPI_Comm comm;
  MPI_Comm_dup(MPI_COMM_WORLD, &comm);

  // run the coupler server and the thermal and flux solvers with different communicators
  if (world_rank < thermal_solver_ranks)
  {
    MPI_Comm_split(comm, 0, world_rank, &thermal_comm);
    int thermal_solver_return_flag = thermal_solver(mesh_file_name, thermal_comm);
    if (thermal_solver_return_flag == 1)
    {
      std::cout << "The thermal solver failed\n";
      // return 1 and abort the program
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Comm_free(&thermal_comm);
  } else if ((world_rank >= thermal_solver_ranks) && 
              (world_rank < thermal_solver_ranks + flux_solver_ranks))
  {
    MPI_Comm_split(comm, 1, world_rank, &flux_comm);
    int flux_solver_return_flag = flux_solver(mesh_file_name, flux_comm);
    if (flux_solver_return_flag == 1)
    {
      std::cout << "The flux solver failed\n";
      // return 1 and abort the program
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Comm_free(&flux_comm);
  } else if ((world_rank < thermal_solver_ranks + flux_solver_ranks + coupler_server_ranks) && 
              (world_rank >= thermal_solver_ranks + flux_solver_ranks))
  {
    ::gmsh::initialize();{
    MPI_Comm_split(comm, 2, world_rank, &server_comm);
    // load the mesh using Omega_h
    auto lib = Omega_h::Library(&argc, &argv, server_comm);
    const auto world = lib.world();
    auto mesh = Omega_h::gmsh::read_parallel("../mesh/parmesh/cylinder", world);
    // hold all the ranks here with MPI_barrier so that the reading is done in parallel and then the
    // coupler server is called
    
    // print the owned entities of the mesh
    //auto nnents_owned = mesh.nents_owned(0);
    auto owned = mesh.owned(0);
    auto nents_owned = std::accumulate(owned.begin(), owned.end(), 0);
    std::stringstream ss;
    ss << "rank: " << world_rank << " The number of owned entities is: " << nents_owned << "\n";
    std::cout << ss.str();

    MPI_Barrier(server_comm);
    mfem_coupler(server_comm, mesh);
    MPI_Comm_free(&server_comm);
    ::gmsh::finalize();}
  } else
  {
    // just do nothing and print that the ranks are not used
    std::cout << "Rank " << world_rank << " not used\n";
  }


  MPI_Comm_free(&comm);
  // finalize MPI
  MPI_Finalize();
  if (world_rank == 0)
  {
    std::cout << "The test passed\n";
  }
  return 0;
}
