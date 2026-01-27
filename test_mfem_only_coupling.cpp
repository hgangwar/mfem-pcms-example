/* ! Description:
* This file implements a verification test for domain decomposition and
 * partitioned coupling strategies using MFEM on a shared 3D mesh.
 *
 * A steady-state heat conduction problem is solved on a cubic domain using
 * linear finite elements. The domain is decomposed into three subregions
 * along the x-direction, and different solve configurations are realized
 * by selectively activating subdomains while preserving global mesh
 * indexing.
 *
 * Three solve modes are supported:
 *   - Client A: subdomains (1,2) active, subdomain 3 inactive
 *   - Client B: subdomains (2,3) active, subdomain 1 inactive
 *   - Client M: all subdomains (1,2,3) active (monolithic reference)
 *
 * Inactive subdomains are numerically decoupled using a near-zero thermal
 * conductivity, allowing partial-domain solves without modifying the mesh.
 * Internal coupling is emulated by imposing Dirichlet constraints on an
 * internal interface plane (x = 0.5 * Lx), where temperature values may be
 * supplied from an external solver or a previous iteration, mimicking
 * data exchange in a coupled multi-physics workflow.
 *
 * Problem setup:
 *   - Governing equation: -div(k(x) grad T) = q(x)
 *   - Heat generation applied only in the central subdomain
 *   - Dirichlet boundary conditions:
 *       T = 300 at x = 0
 *       T = 350 at x = Lx
 *   - All other boundaries are natural (homogeneous Neumann)
 *
 * Numerical details:
 *   - Mesh: structured 3D Cartesian mesh with HEX elements
 *   - Finite elements: H1 Lagrange, order 1
 *   - Linear solver: Conjugate Gradient with Hypre BoomerAMG preconditioner
 *
 * Purpose:
 *   This test case is intended as a controlled demonstration of domain
 *   decomposition, internal interface enforcement, and coupling mechanics
 *   on a shared finite element mesh. It serves as a foundation for future
 *   integration with PCMS-based multi-application coupling and transient
 *   extensions.
 */


#include <iostream>
#include <cmath>
#include <mpi.h>
#include "mfem.hpp"
#include "mfem_field_adapter.h"
#include <Omega_h_mesh.hpp>
#include <pcms/pcms.h>
#include <pcms/types.h>
#include <Omega_h_file.hpp>
#include <Omega_h_for.hpp>
#include <pcms/adapter/omega_h/omega_h_field.h>
#include "test_support.h"
#include "support.h"

using pcms::Copy;
using pcms::GO;
using pcms::Lagrange;
using pcms::make_array_view;
using pcms::MFEMFieldAdapter;
using pcms::OmegaHFieldAdapter;

using namespace std;
namespace ts = test_support;

// Solve one thermal subproblem of the form:
//   (kappa ∇T, ∇v) + beta (T, v) = (rhs_coeff, v)  with essential BCs on
//   ess_bdr
// rhs_coeff is typically beta * T_other_old
double calculate_rms(const std::vector<double> &original_field,
                   const std::vector<double> &updated_field)
{
  const size_t n = original_field.size();
  if (n == 0)
    return 0.0;

  if (updated_field.size() != n)
  {
    std::cerr << "ERROR: Field size mismatch in calculate_rms().\n";
    return 0.0;
  }

  double sum_sq = 0.0;

  for (size_t i = 0; i < n; i++)
  {
    const double diff = updated_field[i] - original_field[i];
    sum_sq += diff * diff;
  }

  double rms = std::sqrt(sum_sq / double(n));

  std::cout << "RMS of (Original field - Updated field): "
            << rms << std::endl;

  return rms;
}

struct Coupling
{
  std::string name;                                  // Coupler name
  std::unique_ptr<pcms::Coupler> cpl;                // Coupler instance
  std::vector<std::string> app_names;                // Attached apps
  std::vector<std::string> field_names;              // Attached fields
  std::map<std::string, pcms::Application*> apps;    // App_name -> pointer
  std::map<std::string, pcms::CoupledField*> fields; // App_field -> pointer
  bool isServer = false;
};

//--------------------------------------------------------------
// Init_Coupler<TAdapter>
//--------------------------------------------------------------
template <typename Adapter_type>
Coupling Init_Coupler(MPI_Comm comm, const std::string& name,
                      const std::vector<std::string>& app_names,
                      const std::vector<std::string>& field_names,
                      bool isServer, const redev::Partition ptn,
                      Adapter_type Adapter)
{
  Coupling cp;
  cp.name = name;
  cp.app_names = app_names;
  cp.field_names = field_names;
  cp.isServer = isServer;

  if (app_names.size() != field_names.size())
    throw std::runtime_error(
      "Mismatch: app_names and field_names must be of the same size.");
  if (isServer)
    cp.cpl = std::make_unique<pcms::Coupler>(name, comm, isServer, ptn);
  else
    cp.cpl = std::make_unique<pcms::Coupler>(name, comm, isServer, ptn);

  for (size_t i = 0; i < app_names.size(); ++i) {
    auto* app = cp.cpl->AddApplication(app_names[i]);
    cp.fields[app_names[i]] = app->AddField(field_names[i], Adapter);
    cp.apps[app_names[i]] = app;
  }
  return cp;
}

static void app_A(MPI_Comm comm, support::ThermalParams params, string solver_type,
                  string prec_type)
{
  int order = 1;
  // Initialize the FEA System
  support::FEMSystem fem =
    support::Init_FEMSystem(comm, params, order, 'A');
  std::string coupler_name = "mfem_coupler";
  std::vector<string> app_name = {"client_A"};
  std::vector<string> field_name = {"temp"};

  // Initialize the MFEM adapter
  auto adapter = MFEMFieldAdapter(app_name[0], *fem.pmesh, *fem.fes, *fem.x);

  // Initialize coupling interface
  auto client =
    Init_Coupler(comm, coupler_name, app_name, field_name, false, {}, adapter);

  // Initialize global comm on the app
  auto gdi = client.apps["client_A"]->Add_GDI<pcms::GO>("global_comm", comm);

  GO flag = 1;
  auto itr = 1;
  do {

    auto curr_field = *fem.x;
    bool use_interior_bc =
      (itr != 1); // No need to apply internal BC for we don't have a soln yet
    auto residual = support::SolveSystem(fem, solver_type, use_interior_bc,
                                              prec_type, 1e-8, 500, 0);
    fem.x->Save("cube_step_1.sol");
    support::SaveParaviewSolution(*fem.mesh, *fem.x, "temperature", "cube_step_1");

    //  Send from A to C
    client.apps["client_A"]->BeginSendPhase();
    client.fields["client_A"]->Send();
    gdi->Send(&residual, "residual", 1);
    client.apps["client_A"]->EndSendPhase();

    // Receive from C to A
    client.apps["client_A"]->BeginReceivePhase();
    client.fields["client_A"]->Receive();
    client.apps["client_A"]->EndReceivePhase();

    // Step sync
    client.apps["client_A"]->BeginReceivePhase();
    auto done = gdi->Receive("done", 1)[0];

    while (!done) {
      sleep(1);
      done = gdi->Receive("flag", 1)[0];
    }
    flag = gdi->Receive("flag", 1)[0];
    client.apps["client_A"]->EndReceivePhase();
    itr++;
    fem.x->Save("cube_step_5.sol");
    support::SaveParaviewSolution(*fem.mesh, *fem.x, "temperature", "cube_step_5");
  } while (false);
}

static void app_B(MPI_Comm comm, support::ThermalParams params , string solver_type,
                  string prec_type)
{
  int order = 1;

  // Initialize the FEA System
  support::FEMSystem fem =
    support::Init_FEMSystem(comm, params, order, 'B');

  std::string coupler_name = "mfem_coupler";
  std::vector<string> app_name = {"client_B"};
  std::vector<string> field_name = {"temp"};

  // Initialize the MFEM adapter
  auto adapter = MFEMFieldAdapter(app_name[0], *fem.pmesh, *fem.fes, *fem.x);
  auto itr = 1;
  auto flag = 1;
  // Initialize coupling interface
  auto client =
    Init_Coupler(comm, coupler_name, app_name, field_name, false, {}, adapter);

  // Initialize global comm on the app
  auto gdi = client.apps["client_B"]->Add_GDI<pcms::GO>("global_comm", comm);
  GO residual = 0;
  do {

    // Receive from C to B
    client.apps["client_B"]->BeginReceivePhase();
    client.fields["client_B"]->Receive();
    client.apps["client_B"]->EndReceivePhase();

    if (itr > 1 && flag == 0)
      break;
    auto residual = support::SolveSystem(fem, solver_type, true, prec_type,
                                              1e-8, 500, 0);
    fem.x->Save("cube_step_3.sol");
    support::SaveParaviewSolution(*fem.mesh, *fem.x, "temperature", "cube_step_3");

    // Send from B to C
    client.apps["client_B"]->BeginSendPhase();
    client.fields["client_B"]->Send();
    gdi->Send(&residual, "residual", 1);
    client.apps["client_B"]->EndSendPhase();

    // Step sync
    client.apps["client_B"]->BeginReceivePhase();
    auto done = gdi->Receive("done", 1)[0];
    while (!done) {
      sleep(1);
      done = gdi->Receive("flag", 1)[0];
    }
    flag = gdi->Receive("flag", 1)[0];
    printf("received flag at B=%d\n", flag);
    client.apps["client_B"]->EndReceivePhase();
    itr++;

  } while (false);
}

void coupler(MPI_Comm comm, support::ThermalParams params)
{
  // Mesh init
  Omega_h::Library lib(nullptr, nullptr, comm);

  int order = 1;
  // Initialize the FEA System
  support::FEMSystem fem =
    support::Init_FEMSystem(comm, params, order, 'M');

  // Define Partition
  redev::LO dim = 3;
  redev::LOs ranks(1);
  std::iota(ranks.begin(), ranks.end(), 0);
  redev::Reals cuts = {0};
  auto partition = redev::Partition{redev::RCBPtn{dim, ranks, cuts}};

  // Coupling labels
  std::string coupler_name = "mfem_coupler";
  std::vector<string> app_names = {"client_A", "client_B"};
  std::vector<string> field_names = {"temp", "temp"};

  // Initialize the MFEM adapter
  auto adapter = MFEMFieldAdapter(coupler_name, *fem.pmesh, *fem.fes, *fem.x);

  // Initialize coupling interface
  auto server =
    Init_Coupler(comm, coupler_name, app_names, field_names, true, partition,
                 adapter);

  // Initialize global comm on the app
  auto gdi_A = server.apps["client_A"]->Add_GDI<pcms::GO>("global_comm", comm);
  auto gdi_B = server.apps["client_B"]->Add_GDI<pcms::GO>("global_comm", comm);

  GO flag = 1; // True to continue
  int itr = 1;
  float tol = 1e-3;
  GO done = 0;
  do {
    // start step
    done = 0;

    // Step start Field
    std::vector<double> field_C(fem.pmesh->GetNV());
    mfem::GridFunction &gf = *fem.x;
    for (int v = 0; v < fem.pmesh->GetNV(); v++)
    {
      field_C[v] = gf(v);
    }

    // Receive from A to C
    server.apps["client_A"]->BeginReceivePhase();
    server.fields["client_A"]->Receive();
    // auto flag = gdi->Receive( "flag", 1)[0];
    auto residual = gdi_A->Receive("residual", 1)[0];
    printf("received residual at coupler from A=%g\n", residual);
    server.apps["client_A"]->EndReceivePhase();

    // --- after update: read new field values
    std::vector<double> field_AC(fem.pmesh->GetNV());
    mfem::GridFunction &gf_AC = *fem.x;
    for (int v = 0; v < fem.pmesh->GetNV(); v++)
    {
      field_AC[v] = gf_AC(v);
    }

    // --- Calculate the rms between Coupler field and A field
    auto rms = calculate_rms(field_C,
                             field_AC); // converting field_C to read<T>

    printf("rms received at coupler:%f\n", rms);
    flag = (rms > tol);
    fem.x->Save("cube_step_2.sol");
    support::SaveParaviewSolution(*fem.mesh, *fem.x, "temperature", "cube_step_2");
    // Send to App B
    server.apps["client_B"]->BeginSendPhase();
    server.fields["client_B"]->Send();
    gdi_B->Send(&flag, "flag", 1);
    gdi_B->Send(&done, "done", 1);
    server.apps["client_B"]->EndSendPhase();

    // Receive from A to C
    server.apps["client_B"]->BeginReceivePhase();
    server.fields["client_B"]->Receive();
    residual = gdi_B->Receive("residual", 1)[0];
    printf("received residual at coupler from B = %g\n", residual);
    server.apps["client_B"]->EndReceivePhase();

    // --- after update: read new field values
    std::vector<double> field_CB(fem.pmesh->GetNV());
    mfem::GridFunction &gf_CB = *fem.x;
    for (int v = 0; v < fem.pmesh->GetNV(); v++)
    {
      field_CB[v] = gf_CB(v);
    }
    // --- Calculate the rms between B field and Coupler field
    rms = calculate_rms(field_AC,
                        field_CB);
    flag = (rms > tol);
    fem.x->Save("cube_step_4.sol");
    support::SaveParaviewSolution(*fem.mesh, *fem.x, "temperature", "cube_step_4");

    server.apps["client_A"]->BeginSendPhase();
    server.fields["client_A"]->Send(); // field send to A
    server.apps["client_A"]->EndSendPhase();

    // Inform A about the step end
    done = 1;
    server.apps["client_A"]->BeginSendPhase();
    gdi_A->Send(&flag, "flag", 1); // Inform A
    gdi_A->Send(&done, "done", 1);
    server.apps["client_A"]->EndSendPhase();
    // Inform B about the step end
    server.apps["client_B"]->BeginSendPhase();
    gdi_B->Send(&flag, "flag", 1);
    gdi_B->Send(&done, "done", 1);
    server.apps["client_B"]->EndSendPhase();

    printf("sent flag %d, with rms %f coupler to A after itr = %d\n", flag, rms,
           itr);
    itr++;
  } while (false);

  std::cout << "The system converged\n";
}

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv); // MPI init
  const auto clientId = atoi(argv[1]);
  REDEV_ALWAYS_ASSERT(clientId >= -1 && clientId <= 1);
  const auto meshFile = argv[2];

  support::ThermalParams params;
  params.size = {1.0, 1.0, 1.0};
  params.ne   = {10, 10, 10};
  params.q_total = 10.0;
  params.kappa   = 1.0;
  params.rho     = 1.0;
  params.cp      = 1.0;
  params.h_flux  = 0.0;
  params.h_conv  = 0.0;
  params.T_conv  = 0.0;
  params.T_dirichlet = 300.0;

  MPI_Comm comm = MPI_COMM_WORLD;
  {
    switch (clientId) {
      case -1: coupler(comm, params); break;
      case 0: app_A(comm, params, argv[3], argv[4]); break;
      case 1: app_B(comm, params, argv[3], argv[4]); break;
      default:
        std::cerr << "Unhandled client id (should be -1, 0,1)\n";
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
  }
  MPI_Finalize();
  return 0;
}
