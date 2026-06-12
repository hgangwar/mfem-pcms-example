#include "include/coupling_support.h"
#include "include/mfem_field_adapter.h"
#include <pcms/pcms.h>
#include <pcms/create_field.h>
#include <pcms/transfer_field2.h>
#include <pcms/utility/types.h>
#include <Omega_h_file.hpp>

using pcms::Copy;
using pcms::GO;
using pcms::Lagrange;
using pcms::make_array_view;
using pcms::MFEMFieldAdapter;
using pcms::OmegaHFieldAdapter;

using namespace std;
typedef pcms::Real dtype;
struct CoordKey
{
  long long ix;
  long long iy;

  bool operator<(const CoordKey& other) const
  {
    if (ix != other.ix) return ix < other.ix;
    return iy < other.iy;
  }
};

static CoordKey MakeCoordKey(double x, double y, double tol)
{
  return {
    static_cast<long long>(std::llround(x / tol)),
    static_cast<long long>(std::llround(y / tol))
  };
}

void CheckCoordinateOverlap(const Omega_h::Mesh& mesh_A,
                            const Omega_h::Mesh& mesh_B,
                            double xmin,
                            double xmax,
                            double ymin,
                            double ymax,
                            double tol)
{
  auto coords_A = Omega_h::HostRead<double>(mesh_A.coords());
  auto coords_B = Omega_h::HostRead<double>(mesh_B.coords());

  std::map<CoordKey, int> A_pts;
  std::map<CoordKey, int> B_pts;

  for (int v = 0; v < mesh_A.nverts(); ++v) {
    const double x = coords_A[2 * v + 0];
    const double y = coords_A[2 * v + 1];

    if (x >= xmin - tol && x <= xmax + tol &&
        y >= ymin - tol && y <= ymax + tol) {
      A_pts[MakeCoordKey(x, y, tol)] = v;
    }
  }

  for (int v = 0; v < mesh_B.nverts(); ++v) {
    const double x = coords_B[2 * v + 0];
    const double y = coords_B[2 * v + 1];

    if (x >= xmin - tol && x <= xmax + tol &&
        y >= ymin - tol && y <= ymax + tol) {
      B_pts[MakeCoordKey(x, y, tol)] = v;
    }
  }

  int common = 0;
  int missing_in_B = 0;
  int missing_in_A = 0;

  std::cout << "\n[COORD_OVERLAP_CHECK]\n";
  std::cout << "  A overlap vertices = " << A_pts.size() << "\n";
  std::cout << "  B overlap vertices = " << B_pts.size() << "\n";

  int shown = 0;

  for (const auto& [key, a_v] : A_pts) {
    auto it = B_pts.find(key);

    if (it != B_pts.end()) {
      ++common;
    } else {
      ++missing_in_B;

      if (shown < 20) {
        std::cout << "  missing in B: A local_v=" << a_v
                  << " x=" << coords_A[2 * a_v + 0]
                  << " y=" << coords_A[2 * a_v + 1]
                  << "\n";
        ++shown;
      }
    }
  }

  shown = 0;

  for (const auto& [key, b_v] : B_pts) {
    auto it = A_pts.find(key);

    if (it == A_pts.end()) {
      ++missing_in_A;

      if (shown < 20) {
        std::cout << "  missing in A: B local_v=" << b_v
                  << " x=" << coords_B[2 * b_v + 0]
                  << " y=" << coords_B[2 * b_v + 1]
                  << "\n";
        ++shown;
      }
    }
  }

  std::cout << "  common coordinates = " << common << "\n";
  std::cout << "  missing in B = " << missing_in_B << "\n";
  std::cout << "  missing in A = " << missing_in_A << "\n";
}
static void app_A(MPI_Comm comm, const std::string mesh_file,
                  support::ThermalParams params, string solver_type,
                  string prec_type)
{
  // Order of fes assumed
  int order = 1;

  // ΩA: [0,0.6]x[0,1]
  mfem::Mesh mesh(mesh_file, 1, 1);
  mfem::ParMesh pmesh(comm, mesh);

  // States
  double T_left = 270.0;
  double left_bdr_x = 1;
  double right_bdr_x = 3;

  // Estimate tolerance for the mesh
  // const double tol= support::DefaultTolX(pmesh);

  // Initialize the FEA System
  support::FEMSystem fem = support::Init_FEMSystem(&pmesh, order, params.kappa);

  int left_bdr_attr = 1;
  int right_bdr_attr = 3;

  mfem::Array<int> ess_bdrA(pmesh.bdr_attributes.Max());
  ess_bdrA = 0;

  ess_bdrA[left_bdr_attr - 1] = 1;
  ess_bdrA[right_bdr_attr - 1] = 1;

  fem.fes->GetEssentialTrueDofs(ess_bdrA, fem.ess_tdofs);

  support::ApplyBoundaryConstantByAttr(pmesh, *fem.x, left_bdr_attr, T_left);
  support::ApplyBoundaryConstantByAttr(pmesh, *fem.x, right_bdr_attr, 280.0);

  // ---- Output ----
  support::OutputPack outA("schwarz_A", pmesh, *fem.fes);

  // Register fields ONCE (pointers must remain valid)
  outA.pvd.RegisterField("T", fem.x);
  outA.pvd.RegisterField("T_exact", &outA.exact);
  outA.pvd.RegisterField("error", &outA.err);

  std::string coupler_name = "mfem_coupler";
  std::vector<string> app_name = {"client_A"};
  std::vector<string> field_name = {"temp"};
  bool use_mask = true;
  pcms::LO attr = 1; // This attribute represents the coupling domain 1 ->
                     // coupling, 2 -> mesh A
  // Initialize the MFEM adapter
  auto adapter =
    MFEMFieldAdapter(app_name[0], *fem.pmesh, *fem.fes, *fem.x, use_mask, attr);

  // Initialize coupling interface
  auto client = support::Init_Coupler(comm, coupler_name, app_name, field_name,
                                      false, {}, adapter);

  // Initialize global comm on the app
  auto gdi = client.apps["client_A"]->Add_GDI<pcms::GO>("global_comm", comm);

  GO flag = 1;
  auto itr = 1;
  do {
    auto curr_field = *fem.x;

    support::PrintTempStats(pmesh, *fem.x, "App A before solve", itr);
    auto residual =
      support::SolveSystem(fem, solver_type, prec_type, 1e-8, 500);
    support::SaveFields(outA, fem, itr);
    support::PrintTempStats(pmesh, *fem.x, "App A after solve: ", itr);

    double err = support::ComputeAbsoluteError(pmesh, *fem.x);

    printf("A abs error=%e\n", err);
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
  } while (flag);
}

static void app_B(MPI_Comm comm, const std::string mesh_file,
                  support::ThermalParams params, string solver_type,
                  string prec_type)
{
  // Order of fes assumed
  int order = 1;
  mfem::Mesh mesh(mesh_file, 1, 1);
  mfem::ParMesh pmesh(comm, mesh);

  // States
  double T_right = 300.0;

  // Estimate tolerance for the mesh
  const double tol = support::DefaultTolX(pmesh);

  // Initialize the FEA System
  support::FEMSystem fem = support::Init_FEMSystem(&pmesh, order, params.kappa);

  // Essential boundaries: both x-min and x-max for each subdomain
  int left_bdr_attr = 1;
  int right_bdr_attr = 3;

  mfem::Array<int> ess_bdrB(pmesh.bdr_attributes.Max());
  ess_bdrB = 0;

  ess_bdrB[left_bdr_attr - 1] = 1;
  ess_bdrB[right_bdr_attr - 1] = 1;

  fem.fes->GetEssentialTrueDofs(ess_bdrB, fem.ess_tdofs);

  // Apply boundary condtition (right edge only)
  support::ApplyBoundaryConstantByAttr(pmesh, *fem.x, right_bdr_attr, T_right);

  // ---- Output ----
  support::OutputPack outB("schwarz_B", pmesh, *fem.fes);

  // Register fields ONCE (pointers must remain valid)
  outB.pvd.RegisterField("T", fem.x);
  outB.pvd.RegisterField("T_exact", &outB.exact);
  outB.pvd.RegisterField("error", &outB.err);

  std::string coupler_name = "mfem_coupler";
  std::vector<string> app_name = {"client_B"};
  std::vector<string> field_name = {"temp"};
  pcms::LO attr = 1; // This indicates the coupling domain
  bool use_mask = true;

  // Initialize the MFEM adapter
  auto adapter =
    MFEMFieldAdapter(app_name[0], *fem.pmesh, *fem.fes, *fem.x, use_mask, attr);
  auto itr = 1;
  auto flag = 1;

  // Initialize coupling interface
  auto client = support::Init_Coupler(comm, coupler_name, app_name, field_name,
                                      false, {}, adapter);

  // Initialize global comm on the app
  auto gdi = client.apps["client_B"]->Add_GDI<pcms::GO>("global_comm", comm);
  GO residual = 0;
  do {
    support::PrintTempStats(pmesh, *fem.x, "App B before receive: ", itr);
    // Receive from C to B
    client.apps["client_B"]->BeginReceivePhase();
    client.fields["client_B"]->Receive();
    client.apps["client_B"]->EndReceivePhase();
    support::ApplyBoundaryConstantByAttr(pmesh, *fem.x, right_bdr_attr,
                                         T_right);

    if (itr > 1 && flag == 0)
      break;
    support::PrintTempStats(pmesh, *fem.x, "App B after receive and before solve: ", itr);

    double err = support::ComputeAbsoluteError(pmesh, *fem.x);

    printf("B before solve abs error=%e\n", err);
    auto residual =
      support::SolveSystem(fem, solver_type, prec_type, 1e-8, 500);
    support::PrintTempStats(pmesh, *fem.x, "App B after solve: ", itr);

    err = support::ComputeAbsoluteError(pmesh, *fem.x);

    printf("B abs error=%e\n", err);
    support::SaveFields(outB, fem, itr);

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
void coupler(MPI_Comm comm, const std::string mesh_A_file,
             const std::string mesh_B_file)
{
  // Mesh init
  Omega_h::Library lib(nullptr, nullptr, comm);
  auto world = lib.world();

  // Read Mesh for App A
  Omega_h::Mesh mesh_A(&lib);
  Omega_h::binary::read(mesh_A_file, world, &mesh_A);

  auto dim = mesh_A.dim();
  const auto nverts = mesh_A.nverts();

  // Read Mesh for App B
  Omega_h::Mesh mesh_B(&lib);
  Omega_h::binary::read(mesh_B_file, world, &mesh_B);

  support::dtype random_temp = 280;
  Omega_h::Read<support::dtype> init(nverts,
                                     random_temp); // init with random guess
  auto field_name = std::string("temp");
  mesh_A.add_tag<support::dtype>(Omega_h::VERT, field_name, 1, init);
  mesh_B.add_tag<support::dtype>(Omega_h::VERT, field_name, 1, init);
  auto isOwned = mesh_A.owned(0);

  pcms::LO tag = 0;
  auto is_overlap_A = support::create_mask(mesh_A, "domain", tag);
  auto is_overlap_B = support::create_mask(mesh_B, "domain", tag);

  // Define Partition
  redev::LOs ranks(1);
  std::iota(ranks.begin(), ranks.end(), 0);
  redev::Reals cuts = {0};
  auto partition = redev::Partition{redev::RCBPtn{dim, ranks, cuts}};

  // Coupling labels
  std::string coupler_name = "mfem_coupler";
  std::vector<string> app_A = {"client_A"};
  std::vector<string> app_B = {"client_B"};
  std::vector<string> field_names = {field_name};

  // Initialize coupling interface
  auto server_A =
    support::Init_Coupler_OH<dtype>(comm, coupler_name, app_A, field_names,
                                    true, partition, mesh_A, is_overlap_A);
  auto server_B =
    support::Init_Coupler_OH<dtype>(comm, coupler_name, app_B, field_names,
                                    true, partition, mesh_B, is_overlap_B);

  // Initialize global comm on the app

  auto gdi_A =
    server_A.apps["client_A"]->Add_GDI<pcms::GO>("global_comm", comm);
  auto gdi_B =
    server_B.apps["client_B"]->Add_GDI<pcms::GO>("global_comm", comm);

  GO flag = 1; // True to continue
  int itr = 1;
  float tol = 1e-3;
  GO done = 0;
  pcms::Real fill_value = 0.0;

  // Setup layouts and  field pointers
  auto layout_A = pcms::CreateLagrangeLayout(
    mesh_A, 1, 1, pcms::CoordinateSystem::Cartesian, "global");
  auto field_A = layout_A->CreateFieldReal();
  field_A->SetOutOfBoundsMode(pcms::OutOfBoundsMode::FILL, fill_value);

  auto layout_B = pcms::CreateLagrangeLayout(
    mesh_B, 1, 1, pcms::CoordinateSystem::Cartesian, "global");
  auto field_B = layout_B->CreateFieldReal();
  field_B->SetOutOfBoundsMode(pcms::OutOfBoundsMode::FILL, fill_value);

  double w = 1.0; // Schwarz coupling relaxation
  do {
    // start step
    done = 0;

    auto dof_C =
      Omega_h::deep_copy(mesh_A.get_array<support::dtype>(0, "temp"));

    // Receive from A to C
    server_A.apps["client_A"]->BeginReceivePhase();
    server_A.fields["client_A"]->Receive();
    auto residual = gdi_A->Receive("residual", 1)[0];
    // printf("received residual at coupler from A=%g\n", residual);
    server_A.apps["client_A"]->EndReceivePhase();

    // Extract the adapter
    const auto& OH_adapter_A = server_A.fields["client_A"]->GetFieldAdapter<pcms::OmegaHFieldAdapter<dtype>>();
    auto& Meshfield_adapter_A = OH_adapter_A->GetField();
    auto nodal_data_A = get_nodal_data(Meshfield_adapter_A);
    for (int i = 0; i < nodal_data_A.size(); i++) {
      printf("Nodal data at index %lu : %f\n", i, nodal_data_A[i]);
    }
    // init the field
    auto dof_A = mesh_A.get_array<support::dtype>(0, "temp");
    double errA = support::ComputeAbsoluteError(mesh_A);
    support::PrintTempStats(mesh_A, "Coupler mesh A", itr);

    printf("iteration=%d A abs error=%e\n", itr, errA);

    // --- after update: read new field values
    auto rms = support::ComputeRMS(Omega_h::Read<support::dtype>(dof_C), dof_A);
    printf("rms received at coupler:%f\n", rms);
    flag = (rms > tol);

    // Interpolate from A to B (source A, target B)
    CheckCoordinateOverlap(mesh_A, mesh_B,
                       0.4, 0.6,
                       0.0, 1.0,
                       1e-10);
    support::PrintTempStats(mesh_B, "Coupler mesh_B before interpolation", itr);
    auto before = Omega_h::deep_copy(mesh_B.get_array<support::dtype>(0, "temp"));
    pcms::interpolate_field2(*field_A, *field_B);
    auto after = mesh_B.get_array<support::dtype>(0, "temp");

    std::cout << "interp rms B = "
              << support::ComputeRMS(before, after)
              << "\n";
    support::PrintTempStats(mesh_B, "Coupler mesh_B after interpolation", itr);

    // Send to App B
    server_B.apps["client_B"]->BeginSendPhase();
    server_B.fields["client_B"]->Send();
    gdi_B->Send(&flag, "flag", 1);
    gdi_B->Send(&done, "done", 1);
    server_B.apps["client_B"]->EndSendPhase();

    // Receive to mesh_B from app B
    server_B.apps["client_B"]->BeginReceivePhase();
    server_B.fields["client_B"]->Receive();
    residual = gdi_B->Receive("residual", 1)[0];
    printf("received residual at coupler from B = %g\n", residual);
    server_B.apps["client_B"]->EndReceivePhase();
    support::PrintTempStats(mesh_B, "Coupler mesh_B after receive back from  app B", itr);

    // Interpolate from A to B (source B, target A)
    support::PrintTempStats(mesh_A, "Coupler mesh_A before interpolation", itr);
    pcms::interpolate_field2(*field_B, *field_A);
    support::PrintTempStats(mesh_A, "Coupler mesh_A after interpolation", itr);

    // --- after update: read new field values
    rms = support::ComputeRMS(Omega_h::Read<support::dtype>(dof_C), dof_A);

    // Send to App A
    server_A.apps["client_A"]->BeginSendPhase();
    server_A.fields["client_A"]->Send();
    gdi_A->Send(&flag, "flag", 1);
    gdi_A->Send(&done, "done", 1);
    server_A.apps["client_A"]->EndSendPhase();

    printf("rms received at coupler:%f\n", rms);
    flag = (rms > tol);

    // Inform A about the step end
    done = 1;
    server_A.apps["client_A"]->BeginSendPhase();
    gdi_A->Send(&flag, "flag", 1); // Inform A
    gdi_A->Send(&done, "done", 1);
    server_A.apps["client_A"]->EndSendPhase();
    // Inform B about the step end
    server_B.apps["client_B"]->BeginSendPhase();
    gdi_B->Send(&flag, "flag", 1);
    gdi_B->Send(&done, "done", 1);
    server_B.apps["client_B"]->EndSendPhase();

    printf("sent flag %d, with rms %f coupler to A after itr = %d\n", flag, rms,
           itr);
    double errB = support::ComputeAbsoluteError(mesh_B);

    printf("iteration=%d B abs error=%e\n", itr, errB);
    support::PrintTempStats(mesh_B, "Coupler mesh B", itr);
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
  params.size = {0.6, 1.0};
  params.ne = {30, 30};
  params.kappa = 1.0;

  MPI_Comm comm = MPI_COMM_WORLD;
  {
    switch (clientId) {
      case -1: coupler(comm, meshFile, argv[3]); break;
      case 0: app_A(comm, meshFile, params, argv[3], argv[4]); break;
      case 1: app_B(comm, meshFile, params, argv[3], argv[4]); break;
    }
  }
  MPI_Finalize();
  return 0;
}
