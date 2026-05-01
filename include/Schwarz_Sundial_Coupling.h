//
// Created by gangwh on 4/22/26.
//

#ifndef PCMS_MFEM_COUPLING_SCHWARTZ_SUNDIAL_COUPLING_H
#define PCMS_MFEM_COUPLING_SCHWARTZ_SUNDIAL_COUPLING_H



#include "mfem.hpp"
#include "support.h"
#include <nvector/nvector_serial.h>
#include <sundials/sundials_stepper.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_errors.h>
#include <sundials/sundials_types.h>

#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace schwarz {

struct Trace
{
  std::vector<double> y;
  std::vector<double> val;

  int Size() const { return static_cast<int>(val.size()); }
  bool Empty() const { return val.empty(); }
};

struct SchwarzConfig
{
  int order = 1;
  int nx = 30;
  int ny = 30;

  double kappa = 1.0;
  double T_left = 270.0;
  double T_right = 300.0;

  double x_A_extract = 0.4;
  double x_A_apply   = 0.6;
  double x_B_apply   = 0.4;
  double x_B_extract = 0.6;

  double omega = 1.0;
  double rel_tol = 1e-12;
  int max_lin_iter = 400;

  double pseudo_dt = 1.0;
  double pseudo_t0 = 0.0;
  double pseudo_tstop = 1e300;

  std::string solver_type = "CG";
  std::string prec_type = "HypreAMG";
};

struct SchwarzStepperContent
{
  SchwarzConfig cfg;

  support::FEMSystem sysA;
  support::FEMSystem sysB;

  int A_attr_xmin = -1;
  int A_attr_xmax = -1;
  int B_attr_xmin = -1;
  int B_attr_xmax = -1;

  double tolA = 1e-12;
  double tolB = 1e-12;

  Trace gA_meta;
  Trace gB_meta;

  sunrealtype tcur = 0.0;
  suncountertype nsteps = 0;
};

void FindXMinMaxBoundaryAttributes(const mfem::ParMesh& pmesh,
                                   int& attr_xmin,
                                   int& attr_xmax);

double DefaultTolX(const mfem::Mesh& mesh);

std::unique_ptr<SchwarzStepperContent>
BuildDefaultSchwarzContent(MPI_Comm comm, const SchwarzConfig& cfg);

// ---------- trace helpers ----------
Trace PairTraceToTrace(const std::vector<std::pair<double, double>>& in);
std::vector<std::pair<double, double>> TraceToPairTrace(const Trace& t);
Trace BlendTrace(const Trace& gA, const Trace& gB, double omega);

std::vector<std::pair<double, double>>
ExtractVertexLineTrace(const mfem::ParMesh& pmesh,
                       const mfem::ParGridFunction& T,
                       double xline,
                       double tol);

void ApplyBoundaryTraceByAttr(mfem::ParMesh& pmesh,
                              mfem::ParGridFunction& gf,
                              int bdr_attr,
                              const std::vector<std::pair<double, double>>& trace,
                              double tol);

void ApplyBoundaryConstantByAttr(mfem::ParMesh& pmesh,
                                 mfem::ParGridFunction& gf,
                                 int bdr_attr,
                                 double value);

void PackState(const Trace& gA, const Trace& gB, N_Vector nv);
void UnpackState(N_Vector nv, Trace& gA, Trace& gB);



// ---------- one Schwarz sweep ----------
int SchwarzSweep(SchwarzStepperContent* C,
                 const Trace& gA_old,
                 const Trace& gB_old,
                 Trace& gA_new,
                 Trace& gB_new);

// ---------- SUNStepper factory ----------
SUNErrCode CreateSchwarzSUNStepper(SUNContext sunctx,
                                   SchwarzStepperContent* content,
                                   SUNStepper* stepper);
SUNErrCode CreateSubdomainASUNStepper(SUNContext sunctx,
                                      SchwarzStepperContent* content,
                                      SUNStepper* stepper);

SUNErrCode CreateSubdomainBSUNStepper(SUNContext sunctx,
                                      SchwarzStepperContent* content,
                                      SUNStepper* stepper);
double RMSDiff(const Trace& a, const Trace& b);

void ErrorToExact_270_30x(const mfem::ParMesh& pmesh,
                          const mfem::ParGridFunction& T,
                          double& rms,
                          double& emax);
} // namespace schwarz

#endif // PCMS_MFEM_COUPLING_SCHWARTZ_SUNDIAL_COUPLING_H