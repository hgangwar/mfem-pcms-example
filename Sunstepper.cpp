//
// Created by gangwh on 4/23/26.
//
#include "Schwarz_Sundial_Coupling.h"

#include <mfem.hpp>
#include <mpi.h>

#include <iostream>
#include <memory>

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  SUNContext sunctx;
  if (SUNContext_Create(SUN_COMM_NULL, &sunctx) != SUN_SUCCESS)
  {
    if (rank == 0) { std::cerr << "SUNContext_Create failed\n"; }
    MPI_Finalize();
    return 1;
  }

  try
  {
    schwarz::SchwarzConfig cfg;
    cfg.omega = 0.7;
    cfg.pseudo_dt = 1.0;
    cfg.pseudo_t0 = 0.0;
    cfg.pseudo_tstop = 100.0;

    auto content_uptr = schwarz::BuildDefaultSchwarzContent(MPI_COMM_WORLD, cfg);
    schwarz::SchwarzStepperContent* content = content_uptr.release();

    const int m = content->gA_meta.Size();
    N_Vector y = N_VNew_Serial(2 * m, sunctx);
    if (!y)
    {
      throw std::runtime_error("Failed to allocate N_Vector state.");
    }

    schwarz::Trace gA0 = content->gA_meta;
    schwarz::Trace gB0 = content->gB_meta;
    for (int i = 0; i < m; i++)
    {
      gA0.val[i] = 282.0;
      gB0.val[i] = 278.0;
    }
    schwarz::PackState(gA0, gB0, y);

    SUNStepper stepper = nullptr;
    if (schwarz::CreateSchwarzSUNStepper(sunctx, content, &stepper) != SUN_SUCCESS)
    {
      N_VDestroy(y);
      delete content;
      throw std::runtime_error("CreateSchwarzSUNStepper failed.");
    }

    schwarz::Trace gA_prev = content->gA_meta;
    schwarz::Trace gB_prev = content->gB_meta;
    schwarz::UnpackState(y, gA_prev, gB_prev);

    const int max_iters = 25;
    const double tol = 1e-8;
    sunrealtype tret = cfg.pseudo_t0;

    for (int k = 0; k < max_iters; k++)
    {
      const sunrealtype tout = tret + cfg.pseudo_dt;
      const SUNErrCode err = SUNStepper_Evolve(stepper, tout, y, &tret);
      if (err != SUN_SUCCESS)
      {
        throw std::runtime_error("SUNStepper_Evolve failed.");
      }

      schwarz::Trace gA = content->gA_meta;
      schwarz::Trace gB = content->gB_meta;
      schwarz::UnpackState(y, gA, gB);

      const double rmsA = RMSDiff(gA, gA_prev);
      const double rmsB = RMSDiff(gB, gB_prev);

      if (rank == 0)
      {
        std::cout << "iter=" << (k + 1)
                  << " tret=" << tret
                  << " rmsA=" << rmsA
                  << " rmsB=" << rmsB << '\n';
      }

      gA_prev = std::move(gA);
      gB_prev = std::move(gB);

      if (std::max(rmsA, rmsB) < tol)
      {
        if (rank == 0)
        {
          std::cout << "Converged: max(rmsA, rmsB) < " << tol << '\n';
        }
        break;
      }
    }

    SUNStepper_Destroy(&stepper);
    N_VDestroy(y);
  }
  catch (const std::exception& e)
  {
    if (rank == 0)
    {
      std::cerr << "ERROR: " << e.what() << '\n';
    }
    SUNContext_Free(&sunctx);
    MPI_Finalize();
    return 1;
  }

  SUNContext_Free(&sunctx);
  MPI_Finalize();
  return 0;
}