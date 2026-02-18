#include "schwartz_coupling_support.h"

using pcms::Copy;
using pcms::GO;
using pcms::Lagrange;
using pcms::make_array_view;
using pcms::MFEMFieldAdapter;
using pcms::OmegaHFieldAdapter;

using namespace std;

OMEGA_H_DEVICE Omega_h::I8 isModelEntInOverlap(const int dim, const int id)
{
  // the TOMMS generated geometric model has
  // entity IDs that increase with the distance
  // from the magnetic axis
  if (dim == 2 && (id >= 22 && id <= 34)) {
    return 1;
  } else if (dim == 1 && (id >= 21 && id <= 34)) {
    return 1;
  } else if (dim == 0 && (id >= 21 && id <= 34)) {
    return 1;
  }
  return 0;
}

/**
 * Create the tag 'isOverlap' for each mesh vertex whose value is 1 if the
 * vertex is classified on a model entity in the closure of the geometric model
 * faces forming the overlap region; the value is 0 otherwise.
 */
Omega_h::Read<Omega_h::I8> markOverlapMeshEntities(Omega_h::Mesh& mesh)
{
  // transfer vtx classification to host
  auto classIds = mesh.get_array<Omega_h::ClassId>(0, "class_id");
  auto classDims = mesh.get_array<Omega_h::I8>(0, "class_dim");
  auto isOverlap = Omega_h::Write<Omega_h::I8>(classIds.size(), "isOverlap");
  auto markOverlap = OMEGA_H_LAMBDA(int i)
  {
    isOverlap[i] = isModelEntInOverlap(classDims[i], classIds[i]);
  };
  Omega_h::parallel_for(classIds.size(), markOverlap);
  auto isOverlap_r = Omega_h::read(isOverlap);
  mesh.add_tag(0, "isOverlap", 1, isOverlap_r);
  return isOverlap_r;
}
Omega_h::HostRead<Omega_h::I8> markMeshOverlapRegion(Omega_h::Mesh& mesh)
{
  auto isOverlap = markOverlapMeshEntities(mesh);
  return Omega_h::HostRead(isOverlap);
}
static void app_A(MPI_Comm comm, support::ThermalParams params, string solver_type,
                  string prec_type)
{
  // Order of fes assumed
  int order = 1;

  // ΩA: [0,0.6]x[0,1]
  mfem::Mesh mesh("/users/gangwh/src/mfem-pcms-example/mesh/box_tri.msh", 1, 1);
  ParMesh pmesh(comm, mesh);

  // States
  double T_left = 270.0;
  double left_bdr_x = 1;
  double right_bdr_x = 3;
  // Estimate tolerance for the mesh
  const double tol= support::DefaultTolX(pmesh);

  // Initialize the FEA System
  support::FEMSystem fem =
    support::Init_FEMSystem(&pmesh, order, params.kappa);

  // Essential boundaries: both x-min and x-max for each subdomain
  Array<int> ess_bdrA(pmesh.bdr_attributes.Max()); ess_bdrA = 0;
  ess_bdrA[left_bdr_x] = 1;
  ess_bdrA[right_bdr_x] = 1;

  // Apply boundary condtition
  fem.fes->GetEssentialTrueDofs(ess_bdrA, fem.ess_tdofs);
  support::ApplyBoundaryConstantByAttr   (pmesh, *fem.x, left_bdr_x, T_left);
  support::ApplyBoundaryConstantByAttr   (pmesh, *fem.x, right_bdr_x, 280);

  // ---- Output ----
  support::OutputPack outA("schwarz_A", pmesh, *fem.fes);

  // Register fields ONCE (pointers must remain valid)
  outA.pvd.RegisterField("T", fem.x);
  outA.pvd.RegisterField("T_exact", &outA.exact);
  outA.pvd.RegisterField("error", &outA.err);

  std::string coupler_name = "mfem_coupler";
  std::vector<string> app_name = {"client_A"};
  std::vector<string> field_name = {"temp"};

  // Initialize the MFEM adapter
  auto adapter = MFEMFieldAdapter(app_name[0], *fem.pmesh, *fem.fes, *fem.x);

  // Initialize coupling interface
  auto client =
    support::Init_Coupler(comm, coupler_name, app_name, field_name, false, {}, adapter);

  // Initialize global comm on the app
  auto gdi = client.apps["client_A"]->Add_GDI<pcms::GO>("global_comm", comm);

  GO flag = 1;
  auto itr = 1;
  do {
    auto curr_field = *fem.x;
    auto residual = support::SolveSystem(fem, solver_type,
                                              prec_type, 1e-8, 500);
    support::SaveFields(outA, fem, itr);
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

static void app_B(MPI_Comm comm, support::ThermalParams params , string solver_type,
                  string prec_type)
{
  // Order of fes assumed
  int order = 1;
  // ΩB: [0.4,1]x[0,1] (build [0,0.6] then shift by +0.4)
  mfem::Mesh mesh("/users/gangwh/src/mfem-pcms-example/mesh/box_tri.msh", 1, 1);
  for (int i = 0; i < mesh.GetNV(); i++) { mesh.GetVertex(i)[0] += 0.4; }
  ParMesh pmesh(comm, mesh);

  // States
  double T_right = 300.0;
  double left_bdr_x = 1;
  double right_bdr_x = 3;
  // Estimate tolerance for the mesh
  const double tol= support::DefaultTolX(pmesh);

  // Initialize the FEA System
  support::FEMSystem fem =
    support::Init_FEMSystem(&pmesh, order, params.kappa);

  // Essential boundaries: both x-min and x-max for each subdomain
  Array<int> ess_bdrB(pmesh.bdr_attributes.Max()); ess_bdrB = 0;
  ess_bdrB[left_bdr_x] = 1;
  ess_bdrB[right_bdr_x] = 1;

  // Set boundary condtition
  fem.fes->GetEssentialTrueDofs(ess_bdrB, fem.ess_tdofs);

  // Apply boundary condtition (left edge should be filled by compiler)
  support::ApplyBoundaryConstantByAttr   (pmesh, *fem.x, left_bdr_x, T_right);

  // ---- Output ----
  support::OutputPack outB("schwarz_B", pmesh, *fem.fes);

  // Register fields ONCE (pointers must remain valid)
  outB.pvd.RegisterField("T", fem.x);
  outB.pvd.RegisterField("T_exact", &outB.exact);
  outB.pvd.RegisterField("error", &outB.err);

  std::string coupler_name = "mfem_coupler";
  std::vector<string> app_name = {"client_B"};
  std::vector<string> field_name = {"temp"};

  // Initialize the MFEM adapter
  auto adapter = MFEMFieldAdapter(app_name[0], *fem.pmesh, *fem.fes, *fem.x);
  auto itr = 1;
  auto flag = 1;

  // Initialize coupling interface
  auto client =
    support::Init_Coupler(comm, coupler_name, app_name, field_name, false, {}, adapter);

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

    auto residual = support::SolveSystem(fem, solver_type, prec_type,
                                              1e-8, 500);

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
void coupler(MPI_Comm comm){
  // Mesh init
  Omega_h::Library lib(nullptr, nullptr, comm);
  auto world = lib.world();
  auto mesh_file = "/users/gangwh/src/mfem-pcms-example/mesh/box_tri.osh";
  // Read Mesh for App A
  Omega_h::Mesh mesh_A(&lib);
  Omega_h::binary::read(mesh_file, world, &mesh_A);

  auto dim = mesh_A.dim();
  const auto nverts = mesh_A.nverts();

  // Create Mesh for App B by shifting mesh_A
  Omega_h::Mesh mesh_B = mesh_A;
  double dx = 0.4;
  support::shift_meshX(mesh_B, dx);

  Omega_h::Write<pcms::Real> init(nverts, 280); // init with random guess
  mesh_A.add_tag<pcms::Real>(Omega_h::VERT, "temp", 1, init);
  mesh_B.add_tag<pcms::Real>(Omega_h::VERT, "temp", 1, init);
  auto isOwned = mesh_A.owned(0);

  // is_overlap is a vector of size mesh.nents(0) and is initialized to 1
  Omega_h::Write<Omega_h::I8> is_overlap(mesh_A.nents(0));
  Omega_h::parallel_for(
    is_overlap.size(), OMEGA_H_LAMBDA(int i) { is_overlap[i] = 1; });

  // Define Partition
  redev::LOs ranks(1);
  std::iota(ranks.begin(), ranks.end(), 0);
  redev::Reals cuts = {0};
  auto partition = redev::Partition{redev::RCBPtn{dim, ranks, cuts}};

  // Coupling labels
  std::string coupler_name = "mfem_coupler";
  std::vector<string> app_names = {"client_A", "client_B"};
  std::vector<string> field_names = {"temp", "temp"};

  // Initialize coupling interface
  auto server_A =
    support::Init_Coupler(comm, coupler_name, app_names, field_names, true, partition,
                 OmegaHFieldAdapter<pcms::Real>("temp", mesh_A, is_overlap));
  auto server_B =
    support::Init_Coupler(comm, coupler_name, app_names, field_names, true, partition,
               OmegaHFieldAdapter<pcms::Real>("temp", mesh_B, is_overlap));

  // Initialize global comm on the app
  auto gdi_A = server_A.apps["client_A"]->Add_GDI<pcms::GO>("global_comm", comm);
  auto gdi_B = server_B.apps["client_B"]->Add_GDI<pcms::GO>("global_comm", comm);

  GO flag = 1; // True to continue
  int itr = 1;
  float tol = 1e-3;
  GO done = 0;

  // Setup layouts and  field pointers
  //auto layout_A = pcms::CreateLagrangeLayout(mesh_A, 1, 1, pcms::CoordinateSystem::Cartesian);
  //auto field_A = layout_A->CreateField();

  //auto layout_B = pcms::CreateLagrangeLayout(mesh_B, 1, 1, pcms::CoordinateSystem::Cartesian);
  //auto field_B = layout_B->CreateField();

  double w = 1.0;  //Schwarz coupling relaxation
  do {
    // start step
    done = 0;

    auto dof_C = Omega_h::deep_copy(mesh_A.get_array<pcms::Real>(0, "temp"));

    // Receive from A to C
    server_A.apps["client_A"]->BeginReceivePhase();
    server_A.fields["client_A"]->Receive();
    auto residual = gdi_A->Receive("residual", 1)[0];
    printf("received residual at coupler from A=%g\n", residual);
    server_A.apps["client_A"]->EndReceivePhase();
    //init the field
    auto dof_A = mesh_A.get_array<pcms::Real>(0, "temp");

    // --- after update: read new field values
    auto rms = support::ComputeRMS(Omega_h::Reals(dof_C),
                             dof_A);
    printf("rms received at coupler:%f\n", rms);
    flag = (rms > tol);

    // Setup gB to be sent to B
    // Extract BC from mesh_A (App A <-> Coupler) to mesh_B (Coupler <-> App_B) gB_new = TB at x=0.6
    auto dof_B = mesh_B.get_array<pcms::Real>(0, "temp");
    auto gB_new = support::ExtractVertexLineTrace(mesh_A, dof_A, 0.6, tol);
    auto gB_exist = support::ExtractVertexLineTrace(mesh_B, dof_B, 0.6, tol);

    // Relax/update gB
    auto gB_relaxed = support::RelaxTrace(gB_exist, gB_new, w);

    // Apply BC to mesh_B by vertex coords
    support::FillTagOnXLineFromTrace(mesh_B, gB_relaxed, 0.6, tol, "temp");

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
    dof_A = mesh_A.get_array<pcms::Real>(0, "temp");

    // Extract BC from mesh_B (Coupler <-> App_B) to mesh_A (App A <-> Coupler)
    dof_B = mesh_B.get_array<pcms::Real>(0, "temp");
    auto gA_new = support::ExtractVertexLineTrace(mesh_B, dof_B, 0.6, tol);
    auto gA_exist = support::ExtractVertexLineTrace(mesh_A, dof_A, 0.6, tol);

    // Relax/update gB
    const auto gA_relaxed = support::RelaxTrace(gA_exist, gA_new, w);

    // --- after update: read new field values
    rms = support::ComputeRMS(Omega_h::Reals(dof_C),
                             dof_A);

    // Apply BC to mesh_A by vertex coords
    support::FillTagOnXLineFromTrace(mesh_A, gA_relaxed, 0.6, tol, "temp");

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
  params.ne   = {30, 30};
  //params.q_total = 10.0;
  params.kappa   = 1.0;
  //params.rho     = 1.0;
  //params.cp      = 1.0;
  //params.h_flux  = 0.0;
  //params.h_conv  = 0.0;
  //params.T_conv  = 0.0;
  //params.T_dirichlet = 300.0;

  MPI_Comm comm = MPI_COMM_WORLD;
  {
    switch (clientId) {
      case -1: coupler(comm); break;
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
