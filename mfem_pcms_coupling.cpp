#include <iostream>
#include <vector>
#include <pcms/pcms.h>
#include <pcms/types.h>
#include "mfem_field_adapter.h"
#include "mfem.hpp"
#include <mpi.h>
#include <Omega_h_mesh.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_for.hpp>
#include <pcms/omega_h_field.h>


using pcms::Copy;
using pcms::CouplerClient;
using pcms::CouplerServer;
using pcms::FieldEvaluationMethod;
using pcms::FieldTransferMethod;
using pcms::GO;
using pcms::Lagrange;
using pcms::make_array_view;
//using pcms::MFEMFieldAdapter;
using pcms::OmegaHFieldAdapter;

//using namespace std::chrono_literals;

//static constexpr bool done = true;
//static constexpr int COMM_ROUNDS = 4;
//namespace ts = test_support;

// This is the coupler of interest
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

  Omega_h::vtk::write_parallel("mfem_couple", &mesh, mesh.dim());

  //std::cout << "mfem coupler function working" << "\n";
}

// main just needs to call mfem_coupler function
int main(int argc, char** argv)

{
  std::cout << "The main function started \n";
  auto lib = Omega_h::Library(&argc, &argv);
  // encapsulates many MPI functions, e.g., world.barrier()
  const auto world = lib.world();
  //create a distributed mesh object, calls mesh.balance() and then returns it
  auto mesh = Omega_h::gmsh::read("../mesh/cylfuelcell3d_full_length.msh", world);
  MPI_Comm comm = world->get_impl();
  
  std::cout << "coupler: mesh is read and coupler is ready to run\n";

  mfem_coupler(comm, mesh);
  std::cout << "here the main program is running \n";
  return 0;
}

