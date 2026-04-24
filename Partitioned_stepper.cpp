//
// Created by gangwh on 4/23/26.
//
// Sunstepper_AB_driver.cpp

#include "Schwarz_Sundial_Coupling.h"

#include <arkode/arkode.h>
#include <arkode/arkode_splittingstep.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_context.h>

#include <mpi.h>
#include <iostream>
#include <stdexcept>
#include <algorithm>

static double ExactTemp(double x)
{
  return 270.0 + 30.0 * x;
}

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  SUNContext sunctx;
  if (SUNContext_Create(SUN_COMM_NULL, &sunctx) != SUN_SUCCESS)
  {
    if (rank == 0) std::cerr << "SUNContext_Create failed\n";
    MPI_Finalize();
    return 1;
  }

  try
  {
    // -----------------------------
    // 1. Build shared Schwarz content
    // -----------------------------
    schwarz::SchwarzConfig cfg;
    cfg.omega = 0.7;
    cfg.pseudo_dt = 1.0;
    cfg.pseudo_t0 = 0.0;
    cfg.pseudo_tstop = 100.0;

    auto content_uptr = schwarz::BuildDefaultSchwarzContent(MPI_COMM_WORLD, cfg);
    auto* content = content_uptr.release();

    const int m = content->gA_meta.Size();

    // State: y = [gA | gB]
    N_Vector y = N_VNew_Serial(2 * m, sunctx);
    if (!y) throw std::runtime_error("Failed to allocate N_Vector.");

    schwarz::Trace gA0 = content->gA_meta;
    schwarz::Trace gB0 = content->gB_meta;

    for (int i = 0; i < m; i++)
    {
      gA0.val[i] = 282.0; // exact T at x=0.4
      gB0.val[i] = 278.0; // initial guess for x=0.6
    }

    schwarz::PackState(gA0, gB0, y);

    // -----------------------------
    // 2. Create two different steppers
    // -----------------------------
    SUNStepper stepperA = nullptr;
    SUNStepper stepperB = nullptr;

    if (schwarz::CreateSubdomainASUNStepper(sunctx, content, &stepperA) != SUN_SUCCESS)
    {
      throw std::runtime_error("CreateSubdomainASUNStepper failed.");
    }

    if (schwarz::CreateSubdomainBSUNStepper(sunctx, content, &stepperB) != SUN_SUCCESS)
    {
      throw std::runtime_error("CreateSubdomainBSUNStepper failed.");
    }

    SUNStepper steppers[2];
    steppers[0] = stepperA;
    steppers[1] = stepperB;

    // -----------------------------
    // 3. Create outer splitting stepper
    // -----------------------------
    void* split_mem = SplittingStepCreate(steppers, 2, cfg.pseudo_t0, y, sunctx);
    if (!split_mem)
    {
      throw std::runtime_error("SplittingStepCreate failed.");
    }

    if (ARKodeSetFixedStep(split_mem, cfg.pseudo_dt) != ARK_SUCCESS)
    {
      throw std::runtime_error("ARKodeSetFixedStep failed.");
    }

    // -----------------------------
    // 4. Evolve pseudo-iterations
    // -----------------------------
    schwarz::Trace gA_prev = content->gA_meta;
    schwarz::Trace gB_prev = content->gB_meta;
    schwarz::UnpackState(y, gA_prev, gB_prev);

    const int max_iters = 25;
    const double tol = 1e-8;

    sunrealtype tret = cfg.pseudo_t0;

    for (int k = 0; k < max_iters; k++)
    {
      const sunrealtype tout = tret + cfg.pseudo_dt;

      const int flag = ARKodeEvolve(split_mem, tout, y, &tret, ARK_NORMAL);
      if (flag < 0)
      {
        throw std::runtime_error("ARKodeEvolve failed.");
      }

      schwarz::Trace gA = content->gA_meta;
      schwarz::Trace gB = content->gB_meta;
      schwarz::UnpackState(y, gA, gB);

      const double rmsA = schwarz::RMSDiff(gA, gA_prev);
      const double rmsB = schwarz::RMSDiff(gB, gB_prev);

      if (rank == 0)
      {
        std::cout << "iter=" << k + 1
                  << " tret=" << tret
                  << " rmsA=" << rmsA
                  << " rmsB=" << rmsB << "\n";
      }

      gA_prev = std::move(gA);
      gB_prev = std::move(gB);

      if (std::max(rmsA, rmsB) < tol)
      {
        if (rank == 0)
        {
          std::cout << "Converged: max trace RMS < " << tol << "\n";
        }
        break;
      }
    }

    // -----------------------------
    // 5. Verify against exact profile
    // -----------------------------
    double rmsEA = 0.0, maxEA = 0.0;
    double rmsEB = 0.0, maxEB = 0.0;

    schwarz::ErrorToExact_270_30x(*content->sysA.pmesh, *content->sysA.x,
                                  rmsEA, maxEA);

    schwarz::ErrorToExact_270_30x(*content->sysB.pmesh, *content->sysB.x,
                                  rmsEB, maxEB);

    if (rank == 0)
    {
      std::cout << "\nFinal solution error against T(x,y)=270+30x\n";
      std::cout << "  Subdomain A: RMS = " << rmsEA
                << "  max = " << maxEA << "\n";
      std::cout << "  Subdomain B: RMS = " << rmsEB
                << "  max = " << maxEB << "\n";
    }

    // -----------------------------
    // 6. Print SUNDIALS stats
    // -----------------------------
    if (rank == 0)
    {
      std::cout << "\nOuter Splitting Stepper Statistics:\n";
      ARKodePrintAllStats(split_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
    }

    // -----------------------------
    // 7. Cleanup
    // -----------------------------
    ARKodeFree(&split_mem);

    // If stepperA/stepperB do not own content separately,
    // avoid double-freeing shared content.
    SUNStepper_Destroy(&stepperA);
    SUNStepper_Destroy(&stepperB);

    N_VDestroy(y);
  }
  catch (const std::exception& e)
  {
    if (rank == 0)
    {
      std::cerr << "ERROR: " << e.what() << "\n";
    }

    SUNContext_Free(&sunctx);
    MPI_Finalize();
    return 1;
  }

  SUNContext_Free(&sunctx);
  MPI_Finalize();

  return 0;
}