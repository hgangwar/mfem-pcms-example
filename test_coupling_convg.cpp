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
// #include <gmsh.h>
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
namespace ts = test_support;
static constexpr bool done = true;
// Solve one thermal subproblem of the form:
//   (kappa ∇T, ∇v) + beta (T, v) = (rhs_coeff, v)  with essential BCs on
//   ess_bdr
// rhs_coeff is typically beta * T_other_old



double calculate_rms(const Omega_h::Read<Omega_h::Real>& original_field,
                     const Omega_h::Read<Omega_h::Real>& updated_field)
{
  const auto n = original_field.size();
  if (n == 0)
    return 0.0;

  Omega_h::Write<Omega_h::Real> sq_diff(n);

  Omega_h::parallel_for(
    n, OMEGA_H_LAMBDA(Omega_h::LO i) {
      const Omega_h::Real diff = updated_field[i] - original_field[i];
      sq_diff[i] = diff * diff;
    });

  const double sum_sq = Omega_h::get_sum(Omega_h::Reals(sq_diff));
  printf("RMS of (Original field - Updated field): %g\n", sum_sq);
  return std::sqrt(sum_sq / static_cast<double>(n));
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

struct FEMSystem
{
  mfem::Mesh* mesh = nullptr;
  mfem::ParMesh* pmesh = nullptr;
  mfem::H1_FECollection* fec = nullptr;
  mfem::ParFiniteElementSpace* fes = nullptr;
  mfem::Array<int> ess_bdr;
  mfem::ParBilinearForm* a = nullptr;
  mfem::ParLinearForm* b = nullptr;
  mfem::ParGridFunction* x = nullptr;
  mfem::DomainLFIntegrator* rhs_int = nullptr;
};

//--------------------------------------------
// 1. Initialization (Mesh + FE + BilinearForm)
//--------------------------------------------
FEMSystem Init_FEMSystem(MPI_Comm comm, const std::string& mesh_file, int order)
{
  FEMSystem sys;

  sys.mesh = new mfem::Mesh(mesh_file.c_str(), 1, 1, true);

  sys.pmesh = new mfem::ParMesh(comm, *sys.mesh);

  int dim = sys.pmesh->Dimension();

  sys.fec = new mfem::H1_FECollection(order, dim);
  sys.fes = new mfem::ParFiniteElementSpace(sys.pmesh, sys.fec);

  sys.ess_bdr.SetSize(sys.pmesh->bdr_attributes.Max());
  sys.ess_bdr = 0;
  if (sys.ess_bdr.Size() > 0)
    sys.ess_bdr[0] = 1;

  sys.a = new mfem::ParBilinearForm(sys.fes);
  mfem::ConstantCoefficient kappa(130.0);
  sys.a->AddDomainIntegrator(new DiffusionIntegrator(kappa));
  sys.a->Assemble();
  sys.a->Finalize();
  mfem::ConstantCoefficient f(100.0);
  sys.b = new mfem::ParLinearForm(sys.fes);
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
                 int max_iter = 500, int print_level = 5)
{
  // --------------------------------------------------
  // Compute essential true DOFs (Dirichlet BCs)
  // --------------------------------------------------
  mfem::Array<int> ess_tdof_list;
  sys.fes->GetEssentialTrueDofs(sys.ess_bdr, ess_tdof_list);

  // --------------------------------------------------
  // Form the parallel linear system A X = B
  // --------------------------------------------------
  mfem::OperatorPtr A; // Owned smart pointer for operator
  mfem::HypreParVector X, B;

  sys.a->FormLinearSystem(ess_tdof_list, *sys.x, *sys.b, A, X, B);

  // Extract actual Hypre matrix pointer for Hypre-based preconditioners
  auto* A_hypre = A.As<mfem::HypreParMatrix>();
  MFEM_VERIFY(A_hypre, "FormLinearSystem did not produce a HypreParMatrix.");

  // --------------------------------------------------
  // Select preconditioner (Hypre-only)
  // --------------------------------------------------
  std::unique_ptr<mfem::Solver> prec;

  if (prec_type == "HypreAMG") {
    auto amg = std::make_unique<mfem::HypreBoomerAMG>(*A_hypre);
    prec = std::move(amg);
  } else if (prec_type == "Jacobi") {
    auto hs = std::make_unique<mfem::HypreSmoother>(*A_hypre);
    hs->SetType(mfem::HypreSmoother::Jacobi);
    prec = std::move(hs);
  } else {
    if (sys.fes->GetParMesh()->GetMyRank() == 0)
      std::cerr << "Unknown preconditioner: " << prec_type << std::endl;
    return;
  }

  // --------------------------------------------------
  // Choose parallel iterative solver
  // --------------------------------------------------
  std::unique_ptr<mfem::IterativeSolver> solver;
  MPI_Comm comm = sys.fes->GetParMesh()->GetComm();

  if (solver_type == "CG")
    solver = std::make_unique<mfem::CGSolver>(comm);
  else if (solver_type == "MINRES")
    solver = std::make_unique<mfem::MINRESSolver>(comm);
  else if (solver_type == "GMRES")
    solver = std::make_unique<mfem::GMRESSolver>(comm);
  else {
    if (sys.fes->GetParMesh()->GetMyRank() == 0)
      std::cerr << "Unknown solver: " << solver_type << std::endl;
    return;
  }

  solver->SetOperator(*A);
  solver->SetPreconditioner(*prec);
  solver->SetRelTol(rel_tol);
  solver->SetAbsTol(0.0);
  solver->SetMaxIter(max_iter);
  solver->SetPrintLevel(print_level);

  // --------------------------------------------------
  // Solve system
  // --------------------------------------------------
  solver->Mult(B, X);

  // --------------------------------------------------
  // Recover finite element solution
  // --------------------------------------------------
  sys.a->RecoverFEMSolution(X, *sys.b, *sys.x);

  if (sys.fes->GetParMesh()->GetMyRank() == 0) {
    std::cout << "Solver converged in " << solver->GetNumIterations()
              << " iterations, final residual = " << solver->GetFinalNorm()
              << std::endl;
  }
}

static void app_A(MPI_Comm comm, string mesh_file, string solver_type,
                  string prec_type)
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

  // Initialize global comm on the app
  auto gdi = client.apps["client_A"]->Add_GDI<pcms::GO>("global_comm", comm);

  GO flag = 1;
  auto itr = 1;
  do {
    auto curr_field = *fem.x;
    SolveSystem(fem, solver_type, prec_type, 1e-8, 500, 0);
    fem.x->Save("cube_step_1.sol");
    auto update_field = *fem.x;
    mfem::ParGridFunction diff(fem.fes);
    diff = update_field.GetData() - curr_field.GetData();
    auto residual = static_cast<int64_t>(diff.Norml2());

    if ( itr >1 && flag == 0)
      break;

    //  Send from A to C
    client.apps["client_A"]->BeginSendPhase();
    client.fields["client_A"]->Send();
    gdi->Send(&residual, "residual", 1);
    printf("Sent flag=%d, residual=%d\n", flag, residual);
    client.apps["client_A"]->EndSendPhase();

    // Receive from C to A
    client.apps["client_A"]->BeginReceivePhase();
    flag = gdi->Receive( "flag", 1)[0];
    printf("received flag=%d, residual=%d\n", flag, residual);
    client.fields["client_A"]->Receive();
    client.apps["client_A"]->EndReceivePhase();

    itr++;

    //fem.x->Save("cube_step_5.sol");
  } while (flag);
}

static void app_B(MPI_Comm comm, string mesh_file, string solver_type,
                  string prec_type)
{
  int order = 1;

  // Initialize the FEA System
  FEMSystem fem = Init_FEMSystem(comm, mesh_file, order);

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

  do {
    fem.x->Save("cube_step_3.sol");
    // Receive from C to B
    client.apps["client_B"]->BeginReceivePhase();
    client.fields["client_B"]->Receive();
    auto flag = gdi->Receive( "flag", 1)[0];
    auto residual = gdi->Receive( "residual", 1)[0];
    client.apps["client_B"]->EndReceivePhase();

    if ( itr > 1 && flag == 0)
      break;
    SolveSystem(fem, solver_type, prec_type, 1e-8, 500, 0);

    //fem.x->Save("cube_step_4.sol");
    // Send from B to C
    client.apps["client_B"]->BeginSendPhase();
    client.fields["client_B"]->Send();
    gdi->Send(&flag, "flag", 1);
    gdi->Send(&residual, "residual", 1);
    client.apps["client_B"]->EndSendPhase();
    itr++;

  } while (flag);
}


void coupler(MPI_Comm comm, std::string mesh_file)
{
  // Mesh init
  Omega_h::Library lib(nullptr, nullptr, comm);
  auto world = lib.world();
  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(mesh_file, world, &mesh);

  // fields init
  const auto nverts = mesh.nverts();
  Omega_h::Write<pcms::Real> init(nverts, 0.0); // init with zero
  mesh.add_tag<pcms::Real>(Omega_h::VERT, "temp", 1, init);
  auto isOwned = mesh.owned(0);

  // is_overlap is a vector of size mesh.nents(0) and is initialized to 1
  Omega_h::Write<Omega_h::I8> is_overlap(mesh.nents(0));
  Omega_h::parallel_for(
    is_overlap.size(), OMEGA_H_LAMBDA(int i) { is_overlap[i] = 1; });

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

  // Initialize coupling interface
  auto server =
    Init_Coupler(comm, coupler_name, app_names, field_names, true, partition,
                 OmegaHFieldAdapter<pcms::Real>("temp", mesh, is_overlap));
  // Initialize global comm on the app
  auto gdi = server.apps["client_A"]->Add_GDI<pcms::GO>("global_comm", comm);
  do {
    auto field_C = Omega_h::deep_copy(mesh.get_array<pcms::Real>(0, "temp"));

    // Receive from A to C
    server.apps["client_A"]->BeginReceivePhase();
    server.fields["client_A"]->Receive();
    auto flag = gdi->Receive( "flag", 1)[0];
    auto residual = gdi->Receive( "residual", 1)[0];
    printf("received flag=%d, residual=%d\n", flag, residual);
    server.apps["client_A"]->EndReceivePhase();

    ts::writeVtk(mesh, "cube_step_", 2);

    // --- after update: read new field values (zero-copy)
    auto field_AC = mesh.get_array<pcms::Real>(0, "temp");

    double rms = calculate_rms(Omega_h::Reals(field_C),
                               field_AC); // converting field_C to read<T>

    // Send to App B
    server.apps["client_B"]->SendPhase(
      [&]() { server.fields["client_B"]->Send(); });

    // Save state before receive from App B
    field_C = Omega_h::deep_copy(field_AC);

    // Receive from A to C
    server.apps["client_B"]->ReceivePhase(
      [&]() { server.fields["client_B"]->Receive(); });

    // --- after update: read new field values (zero-copy)
    auto field_CB = mesh.get_array<pcms::Real>(0, "temp");

    rms = calculate_rms(Omega_h::Reals(field_C),
                        field_CB); // converting field_C to read<T>

    server.apps["client_A"]->BeginSendPhase();
    server.fields["client_A"]->Send();
    gdi->Send(&flag, "flag", 1);
    gdi->Send(&residual, "residual", 1);
    server.apps["client_A"]->EndSendPhase();

  } while (!done);
  std::cout << "The system converged\n";
}

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv); // MPI init
  const auto clientId = atoi(argv[1]);
  REDEV_ALWAYS_ASSERT(clientId >= -1 && clientId <= 1);
  const auto meshFile = argv[2];

  MPI_Comm comm = MPI_COMM_WORLD;
  {
    switch (clientId) {
      case -1: coupler(comm, meshFile); break;
      case 0: app_A(comm, meshFile, argv[3], argv[4]); break;
      case 1: app_B(comm, meshFile, argv[3], argv[4]); break;
      default:
        std::cerr << "Unhandled client id (should be -1, 0,1)\n";
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
  }
  MPI_Finalize();
  return 0;
}
