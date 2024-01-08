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
#define CATCH_CONFIG_RUNNER
#include <iostream>
#include <mpi.h>
//#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include "mfem.hpp"
#include "mfem_field_adapter.h"
#include "mfem_pcms_coupling.cpp"
#include <Omega_h_mesh.hpp>
#include <pcms/pcms.h>
#include <pcms/types.h>
#include <Omega_h_file.hpp>
#include <Omega_h_for.hpp>
#include <redev_variant_tools.h>
#include <pcms/omega_h_field.h>


//using pcms::Copy;
//using pcms::CouplerClient;
//using pcms::CouplerServer;
//using pcms::FieldEvaluationMethod;
//using pcms::FieldTransferMethod;
//using pcms::GO;
//using pcms::Lagrange;
//using pcms::make_array_view;
using pcms::OmegaHFieldAdapter;
using pcms::MFEMFieldAdapter;


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
mfem::GridFunction thermal_solver(const std::string& mesh_file_name, MPI_Comm comm)
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
  std::iota(gf.begin(), gf.end(), 11.23);
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

  // return the grid function
  return gf;
}

// flux solver function for dummy flux solver
// this function is exactly the same as the thermal solver
// except that it receives the data from the thermal solver
// using PCMS
// @param[in] mesh_file_name: the name of the mesh file
// @param[in] comm: MPI communicator
//
// @return the grid function of the flux solver

mfem::GridFunction flux_solver(const std::string& mesh_file_name, MPI_Comm comm)
{
  // load the mesh using mfem
  mfem::Mesh *mesh = new mfem::Mesh(mesh_file_name.c_str(), 1, 1);
  mfem::ParMesh *pmesh = new mfem::ParMesh(comm, *mesh);
  int dim = pmesh->Dimension();

  // create the fe space
  mfem::FiniteElementCollection *fec = new mfem::H1_FECollection(1, dim);
  mfem::ParFiniteElementSpace *fespace = new mfem::ParFiniteElementSpace(pmesh, fec);

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


  // delete the mesh and fe space
  delete fespace;
  delete fec;
  delete pmesh;
  delete mesh;

  // return the grid function
  return gf;
}


// the coupler_run function loads the mesh using Omega_h and
// calls the thermal and flux solvers
// and the coupler server and catches the grid functions
// then compares the data
//

void coupler_run(int argc, char *argv[])
{
  // initialize MPI
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  // mesh name
  std::string mesh_file_name = "../mesh/cylfuelcell3d_full_length.msh";

  // load the mesh using Omega_h
  auto lib = Omega_h::Library(&argc, &argv);
  const auto world = lib.world();
  auto mesh = Omega_h::gmsh::read(mesh_file_name, world);

  // run the couple server
  mfem_coupler(comm, mesh);

  // create the thermal solver
  mfem::GridFunction thermal_gf = thermal_solver(mesh_file_name, comm);

  // create the flux solver
  mfem::GridFunction flux_gf = flux_solver(mesh_file_name, comm);
  // finalize MPI
  MPI_Finalize();
}


// here comes the main function
// for now just check the size  of the grid functions
// TODO: check the data of the grid functions

TEST_CASE("coupling test", "[coupling]")
{
  coupler_run(argc, argv);
  REQUIRE(thermal_gf.Size() == flux_gf.Size());
}
/*
int main(int argc, char *argv[])
{
  // run the test
  int result = Catch::Session().run(argc, argv);
  return result;
}
*/
