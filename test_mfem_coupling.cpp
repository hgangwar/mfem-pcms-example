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


void mfem_coupler(MPI_Comm comm, Omega_h::Mesh& mesh)
{
  // create partition to give it to the server coupler class constructor
  const int dim = 3; // 3D case
  //std::vector<pcms::LO> ranks(8); // does it need to match with # of ranks of comm
  std::vector<pcms::LO> ranks(1); // does it need to match with # of ranks of comm
  std::iota(ranks.begin(),ranks.end(),0);
  //std::vector<pcms::Real> cuts = {0,/*x*/0.5,/*y*/0.75,0.25,/*z*/0.1,0.4,0.8,0.3};
  std::vector<pcms::Real> cuts = {0};
  auto partition = redev::Partition{redev::RCBPtn{dim, ranks, cuts}}; // ? partition constructor work here? How?
  auto& rcb_partition = std::get<redev::RCBPtn>(partition);

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
  //auto* thermal_density_field = thermal_app->AddField(
  //  "density", OmegaHFieldAdapter<GO>("thermal_density", mesh, is_overlap),
  //  FieldTransferMethod::Copy, FieldEvaluationMethod::None,
  //  FieldTransferMethod::Copy, FieldEvaluationMethod::None, is_overlap);
  auto* thermal_density_field = thermal_app->AddField(
    "density", OmegaHFieldAdapter<double>("flux_density", mesh, is_overlap),
    FieldTransferMethod::Copy, FieldEvaluationMethod::None,
    FieldTransferMethod::Copy, FieldEvaluationMethod::None, is_overlap);
  
  std::cout << "coupler: fields are created \n"; 
  // why it is done twice in the proxy_coupling?
  //thermal_app->SendPhase([&]() { thermal_density_field->Send(pcms::Mode::Deferred); });

  std::cout << "coupler: starting receive from thermalSolver\n";
  thermal_app->ReceivePhase([&]() { thermal_density_field->Receive(); });
  std::cout << "coupler: done receiving, starting send\n";
  flux_app->SendPhase([&]() { flux_density_field->Send(); });
  std::cout << "coupler: done sending to fluxSolver\n";

  //Omega_h::vtk::write_parallel("mfem_couple", &mesh, mesh.dim());

  //std::cout << "mfem coupler function working" << "\n";
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
  // fill it with 11.23 value using iota
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
  // fill it with 11.23 value using iota
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
  if (world_size != 3)
  {
    std::cout << "World size must be 3 for " << argv[0] << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // split the communicator based on the rank of each process
  //MPI_Comm comm;
  //int color = world_rank / 3; // Determine color 
  //MPI_Comm_split(MPI_COMM_WORLD, color, world_rank, &comm);


  // mesh name
  std::string mesh_file_name = "../mesh/cylfuelcell3d_full_length.msh";

  // load the mesh using Omega_h
  //auto lib = Omega_h::Library(&argc, &argv, comm);
  //const auto world = lib.world();
  //auto mesh = Omega_h::gmsh::read(mesh_file_name, world);
  
  MPI_Comm server_comm, flux_comm, thermal_comm;
  MPI_Comm comm = MPI_COMM_WORLD;

  // run the coupler server and the thermal and flux solvers with different communicators
  if (world_rank == 0)
  {
    MPI_Comm_split(comm, 0, world_rank, &server_comm);
    // load the mesh using Omega_h
    auto lib = Omega_h::Library(&argc, &argv, server_comm);
    const auto world = lib.world();
    auto mesh = Omega_h::gmsh::read(mesh_file_name, world);
    mfem_coupler(server_comm, mesh);
  } else if (world_rank == 1)
  {
    MPI_Comm_split(comm, 1, world_rank, &flux_comm);
    int flux_solver_return_flag = flux_solver(mesh_file_name, flux_comm);
    if (flux_solver_return_flag == 1)
    {
      std::cout << "The flux solver failed\n";
      // return 1 and abort the program
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  } else if (world_rank == 2)
  {
    MPI_Comm_split(comm, 2, world_rank, &thermal_comm);
    int thermal_solver_return_flag = thermal_solver(mesh_file_name, thermal_comm);
    if (thermal_solver_return_flag == 1)
    {
      std::cout << "The thermal solver failed\n";
      // return 1 and abort the program
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  } 

  // finalize MPI
  MPI_Finalize();
  if (world_rank == 0)
  {
    std::cout << "The test passed\n";
  }
  return 0;
}
