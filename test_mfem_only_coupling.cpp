/* ! Description:
 * This file is a test for the coupling between two MFEM solvers using PCMS.
 * Both of the solvers are based on the same mesh.
 * Loose coupling between TWO steady thermal solvers on a 3D cubic mesh.
 * - Solver A: diffusion with conductivity kappa_A, solved via PCG +
 * Gauss–Seidel
 * - Solver B: diffusion with conductivity kappa_B, solved via MINRES + Jacobi
 * Coupling: symmetric penalty term beta*(T1 - T2) added to both equations;
 * iterate to convergence.
 *
 * Iteration (Gauss–Seidel style):
 *   (kA ∇T1, ∇v) + β (T1, v) = β (T2_old, v)
 *   (kB ∇T2, ∇w) + β (T2, w) = β (T1_new,  w)
 * The coupling is done by PCMS.
 * - Mesh: [0,Lx]x[0,Ly]x[0,Lz] with nx x ny x nz HEX elements
 * - Dirichlet BC: T = T_top on the top face (z = Lz) for both solvers; natural
 * elsewhere.
 * - Coupling parameter beta (>0) controls how strongly T1 and T2 are driven to
 * agree.
 * - Convergence: relative L2 change of T1 and T2 below tolerances.
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
  } while (flag);
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
    fem.x->Save("cube_step_3.sol");
    // Receive from C to B
    client.apps["client_B"]->BeginReceivePhase();
    client.fields["client_B"]->Receive();
    client.apps["client_B"]->EndReceivePhase();

    if (itr > 1 && flag == 0)
      break;
    auto residual = support::SolveSystem(fem, solver_type, true, prec_type,
                                              1e-8, 500, 0);

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

  } while (flag);
}

void coupler(MPI_Comm comm, support::ThermalParams params)
{
  // Mesh init
  Omega_h::Library lib(nullptr, nullptr, comm);
  Mesh serial =
    Mesh::MakeCartesian3D(params.ne[0], params.ne[1], params.ne[2],
                          Element::HEXAHEDRON,
                          params.size[0], params.size[1], params.size[2],
                          true);

  auto mesh  = new Mesh(serial);
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

    //auto field_C = Omega_h::deep_copy(mesh.get_array<pcms::Real>(0, "temp"));
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
    //auto field_AC = mesh.get_array<pcms::Real>(0, "temp");
    std::vector<double> field_AC(fem.pmesh->GetNV());
    mfem::GridFunction &gf_AC = *fem.x;
    for (int v = 0; v < fem.pmesh->GetNV(); v++)
    {
      field_C[v] = gf_AC(v);
    }
    auto rms = calculate_rms(Omega_h::Reals(field_C),
                             field_AC); // converting field_C to read<T>

    printf("rms received at coupler:%f\n", rms);
    flag = (rms > tol);

    // Send to App B
    server.apps["client_B"]->BeginSendPhase();
    server.fields["client_B"]->Send();
    gdi_B->Send(&flag, "flag", 1);
    gdi_B->Send(&done, "done", 1);
    server.apps["client_B"]->EndSendPhase();

    // Save state before receive from App B
    field_C = Omega_h::deep_copy(field_AC);

    // Receive from A to C
    server.apps["client_B"]->BeginReceivePhase();
    server.fields["client_B"]->Receive();
    residual = gdi_B->Receive("residual", 1)[0];
    printf("received residual at coupler from B = %g\n", residual);
    server.apps["client_B"]->EndReceivePhase();

    // --- after update: read new field values
    auto field_CB = mesh.get_array<pcms::Real>(0, "temp");

    rms = calculate_rms(Omega_h::Reals(field_C),
                        field_CB); // converting field_C to read<T>
    flag = (rms > tol);
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
  } while (flag);

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
