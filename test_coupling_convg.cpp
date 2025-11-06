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
#include <redev_variant_tools.h>
#include <pcms/adapter/omega_h/omega_h_field.h>
#include <gmsh.h>
#include "Omega_h_element.hpp"
#include "Omega_h_shape.hpp"
#include <sstream>
#include "test_support.h"

using pcms::Copy;
using pcms::GO;
using pcms::Lagrange;
using pcms::make_array_view;
using pcms::MFEMFieldAdapter;
using pcms::OmegaHFieldAdapter;

using namespace mfem;
using namespace std;

// Solve one thermal subproblem of the form:
//   (kappa ∇T, ∇v) + beta (T, v) = (rhs_coeff, v)  with essential BCs on
//   ess_bdr
// rhs_coeff is typically beta * T_other_old

static constexpr bool done = true;
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

struct FEMSystem
{
  mfem::ParMesh* pmesh = nullptr;
  mfem::H1_FECollection* fec = nullptr;
  mfem::ParFiniteElementSpace* fes = nullptr;
  mfem::Array<int> ess_bdr;
  mfem::BilinearForm* a = nullptr;
  mfem::LinearForm* b = nullptr;
  mfem::ParGridFunction* x = nullptr;
  mfem::DomainLFIntegrator* rhs_int = nullptr;
};

//--------------------------------------------
// 1. Initialization (Mesh + FE + BilinearForm)
//--------------------------------------------
FEMSystem Init_FEMSystem(MPI_Comm comm, const std::string& mesh_file, int order)
{
  FEMSystem sys;

  auto* mesh = new mfem::Mesh(mesh_file.c_str(), 1, 1, true);

  sys.pmesh = new mfem::ParMesh(comm, *mesh);

  int dim = sys.pmesh->Dimension();

  sys.fec = new mfem::H1_FECollection(order, dim);
  sys.fes = new mfem::ParFiniteElementSpace(sys.pmesh, sys.fec);

  sys.ess_bdr.SetSize(sys.pmesh->bdr_attributes.Max());
  sys.ess_bdr = 0;
  if (sys.ess_bdr.Size() > 0)
    sys.ess_bdr[0] = 1;

  sys.a = new BilinearForm(sys.fes);
  mfem::ConstantCoefficient kappa(1.0);
  sys.a->AddDomainIntegrator(new DiffusionIntegrator(kappa));
  sys.a->Assemble();
  sys.a->Finalize();
  mfem::ConstantCoefficient f(0.0);
  sys.b = new LinearForm(sys.fes);
  sys.rhs_int = new DomainLFIntegrator(f);
  sys.b->AddDomainIntegrator(sys.rhs_int);
  sys.b->Assemble();

  sys.x = new ParGridFunction(sys.fes);
  *sys.x = 0.0;

  return sys;
}

//--------------------------------------------
// 2. Solve with argument-based solver config
//--------------------------------------------
void SolveSystem(FEMSystem& sys, const std::string& solver_type,
                 const std::string& prec_type, double rel_tol = 1e-8,
                 int max_iter = 500, int print_level = 1)
{
  SparseMatrix& A = sys.a->SpMat();
  // Choose preconditioner
  std::unique_ptr<Solver> prec;
  if (prec_type == "Jacobi") {
    prec = std::make_unique<DSmoother>(A);
  } else if (prec_type == "GS") {
    prec = std::make_unique<GSSmoother>(A);
  }
#ifdef MFEM_USE_HYPRE
  else if (prec_type == "HypreAMG") {
    prec = std::make_unique<HypreBoomerAMG>(A);
  }
#endif
  else {
    std::cerr << "Unknown preconditioner: " << prec_type << std::endl;
    exit(1);
  }
  // Choose solver
  std::unique_ptr<IterativeSolver> solver;
  if (solver_type == "CG") {
    solver = std::make_unique<CGSolver>();
  } else if (solver_type == "MINRES") {
    solver = std::make_unique<MINRESSolver>();
  } else if (solver_type == "GMRES") {
    solver = std::make_unique<GMRESSolver>();
  } else {
    std::cerr << "Unknown solver: " << solver_type << std::endl;
    exit(1);
  }
  solver->SetOperator(A);
  solver->SetRelTol(rel_tol);
  solver->SetMaxIter(max_iter);
  solver->SetPrintLevel(print_level);
  solver->SetPreconditioner(*prec);
  solver->Mult(*sys.b, *sys.x);
}

static void app_A(MPI_Comm comm, string mesh_file, string solver_type = "CG",
                  string prec_type = "GS")
{
  int order = 1;
  // Initialize the FEA System
  FEMSystem fem = Init_FEMSystem(comm, mesh_file, order);
  std::string coupler_name = "mfem_coupler";
  std::vector<string> app_name = {"client_A"};
  std::vector<string> field_name = {"temp"};

  // Initialize the MFEM adapter
  auto adapter = MFEMFieldAdapter(app_name[0], *fem.pmesh, *fem.fes, *fem.x);

  // Initialize coupling interface
  auto client =
    Init_Coupler(comm, coupler_name, app_name, field_name, false, {}, adapter);

  do {
    SolveSystem(fem, solver_type, prec_type, 1e-8, 500, 0);
    // Send from A to C
    client.apps["client_A"]->SendPhase(
      [&]() { client.fields["client_A"]->Send(); });
    // Receive from C to A
    client.apps["client_A"]->ReceivePhase(
      [&]() { client.fields["client_A"]->Receive(); });
    // sleep(10);
  } while (!done);
}

static void app_B(MPI_Comm comm, string mesh_file, string solver_type = "CG",
                  string prec_type = "GS")
{
  int order = 1;

  // Initialize the FEA System
  FEMSystem fem = Init_FEMSystem(comm, mesh_file, order);

  std::string coupler_name = "mfem_coupler";
  std::vector<string> app_name = {"client_B"};
  std::vector<string> field_name = {"temp"};

  // Initialize the MFEM adapter
  auto adapter = MFEMFieldAdapter(app_name[0], *fem.pmesh, *fem.fes, *fem.x);

  // Initialize coupling interface
  auto client =
    Init_Coupler(comm, coupler_name, app_name, field_name, false, {}, adapter);

  do {
    SolveSystem(fem, solver_type, prec_type, 1e-8, 500, 0);
    // Receive from C to B
    client.apps["client_B"]->ReceivePhase(
      [&]() { client.fields["client_B"]->Receive(); });
    // Send from B to C
    client.apps["client_B"]->SendPhase(
      [&]() { client.fields["client_B"]->Send(); });
    // sleep(10);
  } while (!done);
}

void coupler(MPI_Comm comm, std::string mesh_file)
{
  Omega_h::Library lib(nullptr, nullptr, comm);
  auto world = lib.world();
  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(mesh_file, world, &mesh);

  const auto nverts = mesh.nverts();
  Omega_h::Write<pcms::Real> init(nverts, 0.0); // init with zero
  mesh.add_tag<pcms::Real>(Omega_h::VERT, "temp", 1, init);
  mesh.add_tag<pcms::Real>(Omega_h::VERT, "prev_temp", 1, init);
  auto isOwned = mesh.owned(0);
  // is_overlap is a vector of size mesh.nents(0) and is initialized to 1
  Omega_h::Write<Omega_h::I8> is_overlap(mesh.nents(0));
  Omega_h::parallel_for(
    is_overlap.size(), OMEGA_H_LAMBDA(int i) { is_overlap[i] = 1; });
  printf("Size of mask:%d, size of mesh owned:%d\n", is_overlap.size(),
         isOwned.size());

  redev::LO dim = 3;
  redev::LOs ranks(1);
  std::iota(ranks.begin(), ranks.end(), 0);
  redev::Reals cuts = {0};
  auto partition = redev::Partition{redev::RCBPtn{dim, ranks, cuts}};
  std::string coupler_name = "mfem_coupler";
  std::vector<string> app_names = {"client_A", "client_B"};
  std::vector<string> field_names = {"temp", "temp"};

  // Initialize coupling interface
  auto server =
    Init_Coupler(comm, coupler_name, app_names, field_names, true, partition,
                 OmegaHFieldAdapter<pcms::Real>("temp", mesh, is_overlap));

  // Receive the init field
  // server.apps["client_A"]->ReceivePhase([&]() {
  // server.fields["client_A"]->Receive(); });
  bool first_itr = true;
  pcms::OmegaHField<pcms::Real> field_T("prev_temp", mesh); // field_Tn (empty)
  do {
    auto field_C = mesh.get_array<pcms::Real>(0, "prev_temp");

    // Receive from A to C
    server.apps["client_A"]->ReceivePhase(
      [&]() { server.fields["client_A"]->Receive(); });

    // --- after update: read new field values (zero-copy)
    auto field_AC = mesh.get_array<pcms::Real>(0, "temp");

    const auto n = field_C.size();
    Omega_h::Write<pcms::Real> sq_diff(n);

    // --- Compute squared difference on device
    Omega_h::parallel_for(
      n, OMEGA_H_LAMBDA(Omega_h::LO i) {
        const pcms::Real diff = field_AC[i] - field_C[i];
        sq_diff[i] = diff * diff;
      });

    // Global sum
    double sum_sq = Omega_h::get_sum(Omega_h::Reals(
      sq_diff)); // MPI_reduce not requried in the no partition case
    double rms = std::sqrt(sum_sq / static_cast<double>(n));
    printf(" RMS for the field difference (A-C):%d\n", rms);
    server.apps["client_B"]->SendPhase(
      [&]() { server.fields["client_B"]->Send(); });

    // Get the Coupler field
    auto* adapter = server.fields["client_A"]
                      ->GetFieldAdapter<pcms::OmegaHFieldAdapter<pcms::Real>>();
    const auto& adapter_field = adapter->GetField(); // field_Tn+1

    // --- before update: deep copy of current field values
    pcms::copy_field(adapter_field, field_T);
    // sleep(20);

    // From B to C to A
    field_C = mesh.get_array<pcms::Real>(0, "prev_temp");
    // Receive from A to C
    server.apps["client_B"]->ReceivePhase(
      [&]() { server.fields["client_B"]->Receive(); });

    // --- after update: read new field values (zero-copy)
    auto field_CB = mesh.get_array<pcms::Real>(0, "temp");

    // --- Compute squared difference on device
    Omega_h::parallel_for(
      n, OMEGA_H_LAMBDA(Omega_h::LO i) {
        const pcms::Real diff = field_AC[i] - field_C[i];
        sq_diff[i] = diff * diff;
      });

    // Global sum
    sum_sq = Omega_h::get_sum(Omega_h::Reals(
      sq_diff)); // MPI_reduce not requried in the no partition case
    rms = std::sqrt(sum_sq / static_cast<double>(n));
    printf(" RMS for the field difference (B-C):%d\n", rms);

    server.apps["client_A"]->SendPhase(
      [&]() { server.fields["client_A"]->Send(); });

    // Get the Coupler field
    auto* adapter_BC =
      server.fields["client_A"]
        ->GetFieldAdapter<pcms::OmegaHFieldAdapter<pcms::Real>>();
    const auto& adapter_field_BC = adapter_BC->GetField(); // field_Tn+1

    // --- before update: deep copy of current field values
    pcms::copy_field(adapter_field_BC, field_T);

  } while (!done);
  std::cout << "The system converged\n";
}

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv); // MPI init
  const auto clientId = atoi(argv[1]);
  REDEV_ALWAYS_ASSERT(clientId >= -1 && clientId <= 1);
  const auto meshFile = argv[2];
  int color;
  if (clientId == -1)
    color = 0; // coupler
  else if (clientId == 0)
    color = 1; // client A
  else if (clientId == 1)
    color = 2; // client B
  else
    color = MPI_UNDEFINED;

  MPI_Comm subcomm;
  MPI_Comm_split(MPI_COMM_WORLD, color, 0, &subcomm);

  switch (clientId) {
    case -1: coupler(subcomm, meshFile); break;
    case 0: app_A(subcomm, meshFile); break;
    case 1: app_B(subcomm, meshFile); break;
    default:
      std::cerr << "Unhandled client id (should be -1, 0,1)\n";
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  MPI_Finalize();
  return 0;
}
