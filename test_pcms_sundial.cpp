//
// PCMS outer layer; SUNDIALS only inside coupler
//
// Roles:
//   client_id = -1 : PCMS coupler/server; owns SUNContext, SUNStepper, N_Vector
//   client_id =  0 : App A; MFEM + PCMS only; solve -> send -> receive -> flag
//   client_id =  1 : App B; MFEM + PCMS only; receive -> solve -> send -> flag
//

#include <iostream>
#include <cmath>
#include <numeric>
#include <string>
#include <vector>
#include <unistd.h>
#include <functional>
#include <algorithm>
#include <cstdlib>
#include <map>
#include <memory>
#include <stdexcept>

#include <mpi.h>

#include "mfem.hpp"
#include "include/support.h"
#include "include/mfem_field_adapter.h"
#include "include/test_support.h"
#include "include/Schwarz_Sundial_Coupling.h"

#include <Omega_h_mesh.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_for.hpp>

#include <pcms/pcms.h>
#include <pcms/utility/types.h>

#include <nvector/nvector_serial.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_stepper.h>
#include <sundials/sundials_errors.h>
#include <sundials/sundials_types.h>

typedef pcms::Real dtype;

using pcms::GO;
using pcms::MFEMFieldAdapter;
using pcms::OmegaHFieldAdapter;

namespace ts = test_support;

struct Coupling
{
  std::string name;
  std::unique_ptr<pcms::Coupler> cpl;
  std::vector<std::string> app_names;
  std::vector<std::string> field_names;
  std::map<std::string, pcms::Application*> apps;
  std::map<std::string, pcms::CoupledField*> fields;
  bool isServer = false;
};

template <typename Adapter_type>
Coupling Init_Coupler(MPI_Comm comm,
                      const std::string& name,
                      const std::vector<std::string>& app_names,
                      const std::vector<std::string>& field_names,
                      bool isServer,
                      const redev::Partition ptn,
                      Adapter_type& Adapter)
{
  Coupling cp;
  cp.name = name;
  cp.app_names = app_names;
  cp.field_names = field_names;
  cp.isServer = isServer;

  if (app_names.size() != field_names.size())
  {
    throw std::runtime_error(
        "Mismatch: app_names and field_names must be of the same size.");
  }

  cp.cpl = std::make_unique<pcms::Coupler>(name, comm, isServer, ptn);

  for (size_t i = 0; i < app_names.size(); ++i)
  {
    auto* app = cp.cpl->AddApplication(app_names[i]);
    cp.fields[app_names[i]] = app->AddField(field_names[i], std::move(Adapter));
    cp.apps[app_names[i]] = app;
  }

  return cp;
}

static void CopyNVector(N_Vector src, N_Vector dst)
{
  N_VScale(SUN_RCONST(1.0), src, dst);
}

static double BlockRMSDiff(N_Vector y,
                           N_Vector y_old,
                           int block,
                           int nverts)
{
  const sunrealtype* a = N_VGetArrayPointer(y);
  const sunrealtype* b = N_VGetArrayPointer(y_old);

  const int offset = block * nverts;
  double s = 0.0;

  for (int i = 0; i < nverts; i++)
  {
    const double d = a[offset + i] - b[offset + i];
    s += d * d;
  }

  return std::sqrt(s / std::max(1, nverts));
}

static void PackOmegaHFieldIntoBlock(Omega_h::Mesh& mesh,
                                     const std::string& tag_name,
                                     N_Vector y,
                                     int block,
                                     int nverts)
{
  auto field = mesh.get_array<pcms::Real>(Omega_h::VERT, tag_name);
  Omega_h::HostRead<pcms::Real> hfield(field);

  sunrealtype* data = N_VGetArrayPointer(y);
  const int offset = block * nverts;

  for (int i = 0; i < nverts; i++)
  {
    data[offset + i] = hfield[i];
  }
}

struct CouplerStepperContent
{
  sunrealtype tcur = 0.0;
  sunrealtype dt = 1.0;
  sunrealtype tstop = 1e300;
  suncountertype nsteps = 0;

  double tol = 1e-8;
  bool converged = false;

  int nverts = 0;
  N_Vector yold = nullptr;

  double rmsA = 0.0;
  double rmsB = 0.0;

  std::function<void()> recv_A_field;
  std::function<void()> send_A_field_to_B;
  std::function<void()> recv_B_field;
  std::function<void()> send_B_field_to_A;

  std::function<void(N_Vector)> pack_A_block;
  std::function<void(N_Vector)> pack_B_block;

  std::function<void(GO)> send_flag_to_A;
  std::function<void(GO)> send_flag_to_B;
};

static SUNErrCode CouplerStepper_Evolve(SUNStepper stepper,
                                        sunrealtype tout,
                                        N_Vector y,
                                        sunrealtype* tret)
{
  void* content_void = nullptr;
  if (SUNStepper_GetContent(stepper, &content_void) != SUN_SUCCESS)
  {
    return SUN_ERR_EXT_FAIL;
  }

  auto* C = static_cast<CouplerStepperContent*>(content_void);

  try
  {
    while (C->tcur < tout && C->tcur < C->tstop && !C->converged)
    {
      CopyNVector(y, C->yold);

      // Multiplicative Schwarz / Gauss-Seidel style PCMS exchange:
      // A: solve -> send -> receive -> flag
      // B: receive -> solve -> send -> flag
      C->recv_A_field();
      C->pack_A_block(y);

      C->send_A_field_to_B();

      C->recv_B_field();
      C->pack_B_block(y);

      C->rmsA = BlockRMSDiff(y, C->yold, 0, C->nverts);
      C->rmsB = BlockRMSDiff(y, C->yold, 1, C->nverts);

      const double max_rms = std::max(C->rmsA, C->rmsB);
      C->converged = (max_rms < C->tol);

      GO flag = C->converged ? 0 : 1;

      C->send_B_field_to_A();

      C->send_flag_to_A(flag);
      C->send_flag_to_B(flag);

      C->tcur += C->dt;
      C->nsteps++;
    }

    if (tret) { *tret = C->tcur; }

    SUNStepper_SetLastFlag(stepper, SUN_SUCCESS);
    return SUN_SUCCESS;
  }
  catch (...)
  {
    SUNStepper_SetLastFlag(stepper, SUN_ERR_EXT_FAIL);
    return SUN_ERR_EXT_FAIL;
  }
}

static SUNErrCode CouplerStepper_Reset(SUNStepper stepper,
                                       sunrealtype tR,
                                       N_Vector)
{
  void* content_void = nullptr;
  if (SUNStepper_GetContent(stepper, &content_void) != SUN_SUCCESS)
  {
    return SUN_ERR_EXT_FAIL;
  }

  auto* C = static_cast<CouplerStepperContent*>(content_void);
  C->tcur = tR;
  C->nsteps = 0;
  C->converged = false;
  C->rmsA = 0.0;
  C->rmsB = 0.0;

  return SUN_SUCCESS;
}

static SUNErrCode CouplerStepper_GetNumSteps(SUNStepper stepper,
                                             suncountertype* nsteps)
{
  void* content_void = nullptr;
  if (SUNStepper_GetContent(stepper, &content_void) != SUN_SUCCESS)
  {
    return SUN_ERR_EXT_FAIL;
  }

  auto* C = static_cast<CouplerStepperContent*>(content_void);
  *nsteps = C->nsteps;

  return SUN_SUCCESS;
}

static SUNErrCode CouplerStepper_SetStopTime(SUNStepper stepper,
                                             sunrealtype tstop)
{
  void* content_void = nullptr;
  if (SUNStepper_GetContent(stepper, &content_void) != SUN_SUCCESS)
  {
    return SUN_ERR_EXT_FAIL;
  }

  auto* C = static_cast<CouplerStepperContent*>(content_void);
  C->tstop = tstop;

  return SUN_SUCCESS;
}

static SUNErrCode CouplerStepper_Destroy(SUNStepper stepper)
{
  void* content_void = nullptr;
  if (SUNStepper_GetContent(stepper, &content_void) != SUN_SUCCESS)
  {
    return SUN_ERR_EXT_FAIL;
  }

  delete static_cast<CouplerStepperContent*>(content_void);
  SUNStepper_SetContent(stepper, nullptr);

  return SUN_SUCCESS;
}

static SUNErrCode CreateCouplerSUNStepper(SUNContext sunctx,
                                          CouplerStepperContent* content,
                                          SUNStepper* stepper)
{
  SUNErrCode err = SUNStepper_Create(sunctx, stepper);
  if (err != SUN_SUCCESS) { return err; }

  err = SUNStepper_SetContent(*stepper, content);
  if (err != SUN_SUCCESS) { return err; }

  err = SUNStepper_SetEvolveFn(*stepper, CouplerStepper_Evolve);
  if (err != SUN_SUCCESS) { return err; }

  err = SUNStepper_SetResetFn(*stepper, CouplerStepper_Reset);
  if (err != SUN_SUCCESS) { return err; }

  err = SUNStepper_SetReInitFn(*stepper, CouplerStepper_Reset);
  if (err != SUN_SUCCESS) { return err; }

  err = SUNStepper_SetStopTimeFn(*stepper, CouplerStepper_SetStopTime);
  if (err != SUN_SUCCESS) { return err; }

  err = SUNStepper_SetGetNumStepsFn(*stepper, CouplerStepper_GetNumSteps);
  if (err != SUN_SUCCESS) { return err; }

  err = SUNStepper_SetDestroyFn(*stepper, CouplerStepper_Destroy);
  if (err != SUN_SUCCESS) { return err; }

  return SUN_SUCCESS;
}

static void app_A(MPI_Comm comm,
                  const std::string& mesh_file,
                  const std::string& solver_type,
                  const std::string& prec_type)
{
  constexpr int order = 1;
  constexpr double kappa = 1.0;

  mfem::Mesh mesh(mesh_file.c_str(), 1, 1);
  mfem::ParMesh pmesh(comm, mesh);

  auto fem = support::Init_FEMSystem(&pmesh, order, kappa);

  const std::string coupler_name = "mfem_coupler";
  const std::vector<std::string> app_names = {"client_A"};
  const std::vector<std::string> field_names = {"temp"};

  auto adapter = MFEMFieldAdapter("client_A", *fem.pmesh, *fem.fes, *fem.x);
  auto client =
      Init_Coupler(comm, coupler_name, app_names, field_names, false, {}, adapter);

  auto gdi = client.apps["client_A"]->Add_GDI<pcms::GO>("global_comm", comm);

  GO flag = 1;
  int itr = 1;

  while (flag)
  {
    // A sends first.
    const auto residual =
        support::SolveSystem(fem, solver_type, prec_type, 1e-8, 500);

    client.apps["client_A"]->BeginSendPhase();
    client.fields["client_A"]->Send();
    client.apps["client_A"]->EndSendPhase();

    // Then A receives B-updated field from the coupler.
    client.apps["client_A"]->BeginReceivePhase();
    client.fields["client_A"]->Receive();
    client.apps["client_A"]->EndReceivePhase();

    // Then receive continue/stop flag.
    client.apps["client_A"]->BeginReceivePhase();
    flag = gdi->Receive("flag", 1)[0];
    client.apps["client_A"]->EndReceivePhase();

    std::cout << "[client_A] itr=" << itr
              << " residual=" << residual
              << " flag=" << flag << "\n";

    itr++;
  }

  support::DestroyFEMSystem(fem);
}

static void app_B(MPI_Comm comm,
                  const std::string& mesh_file,
                  const std::string& solver_type,
                  const std::string& prec_type)
{
  constexpr int order = 1;
  constexpr double kappa = 1.0;

  mfem::Mesh mesh(mesh_file.c_str(), 1, 1);
  mfem::ParMesh pmesh(comm, mesh);

  auto fem = support::Init_FEMSystem(&pmesh, order, kappa);

  const std::string coupler_name = "mfem_coupler";
  const std::vector<std::string> app_names = {"client_B"};
  const std::vector<std::string> field_names = {"temp"};

  auto adapter = MFEMFieldAdapter("client_B", *fem.pmesh, *fem.fes, *fem.x);
  auto client =
      Init_Coupler(comm, coupler_name, app_names, field_names, false, {}, adapter);

  auto gdi = client.apps["client_B"]->Add_GDI<pcms::GO>("global_comm", comm);

  GO flag = 1;
  int itr = 1;

  while (flag)
  {
    // B receives first.
    client.apps["client_B"]->BeginReceivePhase();
    client.fields["client_B"]->Receive();
    client.apps["client_B"]->EndReceivePhase();

    // Then B solves and sends.
    const auto residual =
        support::SolveSystem(fem, solver_type, prec_type, 1e-8, 500);

    client.apps["client_B"]->BeginSendPhase();
    client.fields["client_B"]->Send();
    client.apps["client_B"]->EndSendPhase();

    // Then receive continue/stop flag.
    client.apps["client_B"]->BeginReceivePhase();
    flag = gdi->Receive("flag", 1)[0];
    client.apps["client_B"]->EndReceivePhase();

    std::cout << "[client_B] itr=" << itr
              << " residual=" << residual
              << " flag=" << flag << "\n";

    itr++;
  }

  support::DestroyFEMSystem(fem);
}

static void coupler(MPI_Comm comm,
                    const std::string& mesh_file,
                    double tol = 1e-8,
                    int max_iters = 200)
{
  Omega_h::Library lib(nullptr, nullptr, comm);
  auto world = lib.world();

  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(mesh_file, world, &mesh);

  const int nverts = mesh.nverts();

  Omega_h::Write<pcms::Real> init(nverts, 0.0);
  mesh.add_tag<pcms::Real>(Omega_h::VERT, "temp", 1, init);

  Omega_h::Write<Omega_h::I8> is_overlap(mesh.nents(Omega_h::VERT));
  Omega_h::parallel_for(
      is_overlap.size(),
      OMEGA_H_LAMBDA(int i) { is_overlap[i] = 1; });

  redev::LO dim = 3;
  redev::LOs ranks(1);
  std::iota(ranks.begin(), ranks.end(), 0);
  redev::Reals cuts = {0};
  auto partition = redev::Partition{redev::RCBPtn{dim, ranks, cuts}};

  const std::string coupler_name = "mfem_coupler";
  const std::vector<std::string> app_names = {"client_A", "client_B"};
  const std::vector<std::string> field_names = {"temp", "temp"};

  auto adapter = OmegaHFieldAdapter<dtype>("temp", mesh, is_overlap);
  auto server =
      Init_Coupler(comm, coupler_name, app_names, field_names, true, partition, adapter);

  auto gdi_A = server.apps["client_A"]->Add_GDI<pcms::GO>("global_comm", comm);
  auto gdi_B = server.apps["client_B"]->Add_GDI<pcms::GO>("global_comm", comm);

  SUNContext sunctx;
  if (SUNContext_Create(SUN_COMM_NULL, &sunctx) != SUN_SUCCESS)
  {
    throw std::runtime_error("SUNContext_Create failed in coupler.");
  }

  N_Vector y = N_VNew_Serial(2 * nverts, sunctx);
  N_Vector y_old = N_VClone(y);
  if (!y || !y_old)
  {
    throw std::runtime_error("Failed to allocate coupler N_Vector state.");
  }

  N_VConst(SUN_RCONST(0.0), y);
  N_VConst(SUN_RCONST(0.0), y_old);

  auto* stepper_content = new CouplerStepperContent;
  stepper_content->dt = 1.0;
  stepper_content->tcur = 0.0;
  stepper_content->tstop = static_cast<sunrealtype>(max_iters);
  stepper_content->tol = tol;
  stepper_content->nverts = nverts;
  stepper_content->yold = y_old;

  stepper_content->recv_A_field = [&]() {
    server.apps["client_A"]->BeginReceivePhase();
    server.fields["client_A"]->Receive();
    server.apps["client_A"]->EndReceivePhase();
  };

  stepper_content->send_A_field_to_B = [&]() {
    server.apps["client_B"]->BeginSendPhase();
    server.fields["client_B"]->Send();
    server.apps["client_B"]->EndSendPhase();
  };

  stepper_content->recv_B_field = [&]() {
    server.apps["client_B"]->BeginReceivePhase();
    server.fields["client_B"]->Receive();
    server.apps["client_B"]->EndReceivePhase();
  };

  stepper_content->send_B_field_to_A = [&]() {
    server.apps["client_A"]->BeginSendPhase();
    server.fields["client_A"]->Send();
    server.apps["client_A"]->EndSendPhase();
  };

  stepper_content->pack_A_block = [&](N_Vector yy) {
    PackOmegaHFieldIntoBlock(mesh, "temp", yy, 0, nverts);
  };

  stepper_content->pack_B_block = [&](N_Vector yy) {
    PackOmegaHFieldIntoBlock(mesh, "temp", yy, 1, nverts);
  };

  stepper_content->send_flag_to_A = [&](GO flag) {
    server.apps["client_A"]->BeginSendPhase();
    gdi_A->Send(&flag, "flag", 1);
    server.apps["client_A"]->EndSendPhase();
  };

  stepper_content->send_flag_to_B = [&](GO flag) {
    server.apps["client_B"]->BeginSendPhase();
    gdi_B->Send(&flag, "flag", 1);
    server.apps["client_B"]->EndSendPhase();
  };

  SUNStepper stepper = nullptr;
  if (CreateCouplerSUNStepper(sunctx, stepper_content, &stepper) != SUN_SUCCESS)
  {
    throw std::runtime_error("CreateCouplerSUNStepper failed.");
  }

  int itr = 1;
  sunrealtype tret = 0.0;

  while (!stepper_content->converged && itr <= max_iters)
  {
    const sunrealtype tout = tret + stepper_content->dt;

    const SUNErrCode err = SUNStepper_Evolve(stepper, tout, y, &tret);
    if (err != SUN_SUCCESS)
    {
      throw std::runtime_error("Coupler SUNStepper_Evolve failed.");
    }

    std::cout << "[coupler] itr=" << itr
              << " tret=" << tret
              << " rms_A_sun=" << stepper_content->rmsA
              << " rms_B_sun=" << stepper_content->rmsB
              << " max_rms_sun="
              << std::max(stepper_content->rmsA, stepper_content->rmsB)
              << " converged=" << stepper_content->converged
              << "\n";

    itr++;
  }

  if (!stepper_content->converged)
  {
    std::cout << "[coupler] stopped at max_iters=" << max_iters << "\n";
  }
  else
  {
    std::cout << "[coupler] converged based on SUNDIALS DOF state\n";
  }

  SUNStepper_Destroy(&stepper);
  N_VDestroy(y_old);
  N_VDestroy(y);
  SUNContext_Free(&sunctx);
}

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);

  if (argc < 3)
  {
    std::cerr << "Usage:\n"
              << "  " << argv[0] << " -1 mesh.osh\n"
              << "  " << argv[0] << "  0 mesh.mesh solver prec\n"
              << "  " << argv[0] << "  1 mesh.mesh solver prec\n";
    MPI_Finalize();
    return 1;
  }

  const int client_id = std::atoi(argv[1]);
  const std::string mesh_file = argv[2];

  try
  {
    MPI_Comm comm = MPI_COMM_WORLD;

    if (client_id == -1)
    {
      coupler(comm, mesh_file);
    }
    else if (client_id == 0)
    {
      if (argc < 5)
      {
        throw std::runtime_error("client_A needs mesh solver prec.");
      }
      app_A(comm, mesh_file, argv[3], argv[4]);
    }
    else if (client_id == 1)
    {
      if (argc < 5)
      {
        throw std::runtime_error("client_B needs mesh solver prec.");
      }
      app_B(comm, mesh_file, argv[3], argv[4]);
    }
    else
    {
      throw std::runtime_error("client_id must be -1, 0, or 1.");
    }
  }
  catch (const std::exception& e)
  {
    std::cerr << "ERROR: " << e.what() << "\n";
    MPI_Finalize();
    return 1;
  }

  MPI_Finalize();
  return 0;
}