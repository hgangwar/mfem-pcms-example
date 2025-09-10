/* 
This solver solves the neutron diffusion equation in a full length fuel rod.
The fuel rod is a cylinder with radius 0.74 cm and height 140 cm.
Clad and coolant are not included in this model.
*/

// * This is a modified version of the example 11p code from MFEM

/* Description:  Here in this solver we are trying to solve the generalized
               eigenvalue problem -Delta u = B * lambda u with homogeneous
               Dirichlet boundary conditions.

               We only calculate the first eigenvalue and eigenvector.
*/

/** 
* The neutron diffusion equation is given by:
* * -D Delta phi + (Sigma_a - nu Sigma_f) phi = 0
* where:
* * D is the diffusion coefficient
* * Sigma_a is the absorption cross section
* * Sigma_f is the fission cross section
* * nu is the number of neutrons produced per fission
* * phi is the neutron flux

* The neutron diffusion equation can be modified to:
* * Delta phi - lambda B phi = 0
* where:
* * lambda = c2/k - c1 => k = c2/(lambda + c1)
* * B = rho^2 [rho is normalized fuel density: rho = rho_fuel/rho_nominal]
* * c1 = (3 Sigma_a Sigma_t^2) / sigma_s
* * c2 = (3 nu Sigma_f Sigma_t^2) / sigma_s
*/
#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include "pcms/pcms.h"
#include "mfem_field_adapter.h"

using pcms::Coupler;
using pcms::MFEMFieldAdapter;
//# define RAND_MAX 100

using namespace std;
using namespace mfem;

/* ! it's implemented in the mfem adaptedkKDKDkddkr
std::vector<long int> getGids(const ParMesh& pmesh, const ParFiniteElementSpace& fes )
{
	auto * R = fes.GetRestrictionMatrix();
	if(!R) {
          std::cerr<<"R matrix is nullptr\n";
          std::abort();
	}
	Array<HYPRE_BigInt> gids;
	pmesh.GetGlobalVertexIndices(gids);
	int size = gids.Size();
        mfem::Vector gid_vector(size);
        for(int i=0; i<size; ++i) {
          gid_vector[i] = gids[i];
        }
	mfem::Vector tgids(fes.GetTrueVSize());
	//R->BooleanMult(gids, tgids);
	R->Mult(gid_vector, tgids);
	//auto gids_host = gids.HostRead();
	return {tgids.begin(), tgids.end()};
}
*/


int main(int argc, char *argv[])
{
   // 1. Initialize MPI and HYPRE.
   Mpi::Init(argc, argv);
   int num_procs = Mpi::WorldSize();
   int myid = Mpi::WorldRank();
   Hypre::Init();


   //*  Neutronic parameter and calculations
   /** 
   * This parameters are for nominal fuel density
   * @param D Diffusion coefficient
   * @param nu_sigma_f Fission cross section
   * @param sigma_t Total cross section
   * @param sigma_s Scattering cross section
   * @param rho Fuel density (nominal)
   */
   /**
    * ! The one group constants are found to be(using OpenMC):
  	Name	      Mean	   SD
   Total	      0.257243	0.000042
   Absorption	0.013829	0.000002
   Scattering	0.243414	0.000040
   nuFission	0.034901	0.000007
   */

   //double D = 1.3;
   double sigma_a = 0.013829;
   double nu_sigma_f = 0.034901;
   double sigma_t = 0.257243;
   double sigma_s = 0.243414;
   double rho = 11880;

   double c1 = (3 * sigma_t * sigma_t * sigma_a) / sigma_s;
   double c2 = (3 * nu_sigma_f * sigma_t * sigma_t) / sigma_s;


   // 2. Parse command-line options.

   // * We will use this defauld mesh file
   const char *mesh_file = "../mesh/cylfuelcell3d_full_length.msh";

   //const char *mesh_file = "../mesh/testmesh.msh";

   int ser_ref_levels = 0;
   int par_ref_levels = 0;
   int order = 1;
   int nev = 1;
   int seed = 75;
   bool slu_solver  = false;
   bool sp_solver = false;
   bool cpardiso_solver = false;
   bool visualization = 0;
   bool paraview = true;

   OptionsParser args(argc, argv);
   
   // TODO: check for if mesh file is given
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&ser_ref_levels, "-rs", "--refine-serial",
                  "Number of times to refine the mesh uniformly in serial.");
   args.AddOption(&par_ref_levels, "-rp", "--refine-parallel",
                  "Number of times to refine the mesh uniformly in parallel.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree) or -1 for"
                  " isoparametric space.");
   args.AddOption(&nev, "-n", "--num-eigs",
                  "Number of desired eigenmodes.");
   args.AddOption(&seed, "-s", "--seed",
                  "Random seed used to initialize LOBPCG.");


#ifdef MFEM_USE_SUPERLU
   args.AddOption(&slu_solver, "-slu", "--superlu", "-no-slu",
                  "--no-superlu", "Use the SuperLU Solver.");
#endif
#ifdef MFEM_USE_STRUMPACK
   args.AddOption(&sp_solver, "-sp", "--strumpack", "-no-sp",
                  "--no-strumpack", "Use the STRUMPACK Solver.");
#endif
#ifdef MFEM_USE_MKL_CPARDISO
   args.AddOption(&cpardiso_solver, "-cpardiso", "--cpardiso", "-no-cpardiso",
                  "--no-cpardiso", "Use the MKL CPardiso Solver.");
#endif
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.Parse();
   if (slu_solver && sp_solver)
   {
      if (myid == 0)
         cout << "WARNING: Both SuperLU and STRUMPACK have been selected,"
              << " please choose either one." << endl
              << "         Defaulting to SuperLU." << endl;
      sp_solver = false;
   }
   // The command line options are also passed to the STRUMPACK
   // solver. So do not exit if some options are not recognized.
   if (!sp_solver)
   {
      if (!args.Good())
      {
         if (myid == 0)
         {
            args.PrintUsage(cout);
         }
         return 1;
      }
   }
   if (myid == 0)
   {
      args.PrintOptions(cout);
   }

   // 3. Read the (serial) mesh from the given mesh file on all processors. We
   //    can handle triangular, quadrilateral, tetrahedral, hexahedral, surface
   //    and volume meshes with the same code.
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   int dim = mesh->Dimension();

   // 4. Refine the serial mesh on all processors to increase the resolution. In
   //    this example we do 'ref_levels' of uniform refinement (2 by default, or
   //    specified on the command line with -rs).
   for (int lev = 0; lev < ser_ref_levels; lev++)
   {
      mesh->UniformRefinement();
   }

   // 5. Define a parallel mesh by a partitioning of the serial mesh. Refine
   //    this mesh further in parallel to increase the resolution (1 time by
   //    default, or specified on the command line with -rp). Once the parallel
   //    mesh is defined, the serial mesh can be deleted.
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
   delete mesh;
   for (int lev = 0; lev < par_ref_levels; lev++)
   {
      pmesh->UniformRefinement();
   }

   // 6. Define a parallel finite element space on the parallel mesh. Here we
   //    use continuous Lagrange finite elements of the specified order. If
   //    order < 1, we instead use an isoparametric/isogeometric space.
   FiniteElementCollection *fec;
   if (order > 0)
   {
      fec = new H1_FECollection(order, dim);
   }
   else if (pmesh->GetNodes())
   {
      fec = pmesh->GetNodes()->OwnFEC();
   }
   else
   {
      fec = new H1_FECollection(order = 1, dim);
   }
   ParFiniteElementSpace *fespace = new ParFiniteElementSpace(pmesh, fec);
   HYPRE_BigInt size = fespace->GlobalTrueVSize();
   if (myid == 0)
   {
      cout << "Number of unknowns: " << size << endl;
   }


   // receive density from the thermal solver
   ParGridFunction dent(fespace);
   Coupler cpl("fluxClient", MPI_COMM_WORLD);
   // ? No addfield needed? how it's gonna know where to save the data?
   cpl.AddField("density", MFEMFieldAdapter(std::string("flux2th"), *pmesh, *fespace, dent));

   cpl.BeginReceivePhase();
   cpl.ReceiveField("temp");
   cpl.EndReceivePhase();


  /* 
   // Call the getGid function and print
   auto gids = getGids(*pmesh, *fespace);
   std::stringstream ss;
   ss << myid<<": ";
   for (auto id : gids)
   {
           ss << id << " ";
   }
   ss<<"\n";
   std::cout<<ss.str();
  */

   // 7. Set up the parallel bilinear forms a(.,.) and m(.,.) on the finite
   //    element space. The first corresponds to the Laplacian operator -Delta,
   //    while the second is a simple mass matrix needed on the right hand side
   //    of the generalized eigenvalue problem below. The boundary conditions
   //    are implemented by elimination with special values on the diagonal to
   //    shift the Dirichlet eigenvalues out of the computational range. After
   //    serial and parallel assembly we extract the corresponding parallel
   //    matrices A and M.
   ConstantCoefficient one(1.0);
   // ! FunctionCoefficient ff(myF);


   Array<int> ess_bdr;
   if (pmesh->bdr_attributes.Size())
   {
      ess_bdr.SetSize(pmesh->bdr_attributes.Max());
      ess_bdr = 1; // boundary condition for the fuel rod
   }

   ParBilinearForm *a = new ParBilinearForm(fespace);
   a->AddDomainIntegrator(new DiffusionIntegrator(one));
   if (pmesh->bdr_attributes.Size() == 0)
   {
      // Add a mass term if the mesh has no boundary, e.g. periodic mesh or
      // closed surface.
      a->AddDomainIntegrator(new MassIntegrator(one));
   }
   a->Assemble();
   a->EliminateEssentialBCDiag(ess_bdr, 1.0);
   a->Finalize();


   // * Bilinear form for the mass matrix
   // we need density for the mass matrix
   // this grid function will come from the temperature solver
   // here we will normalize the density to the nominal density and square it to get the B matrix
   ParGridFunction density(fespace);
   // this will be replaced when we get the density from the temperature solver
   ConstantCoefficient rho_coeff(rho); // rho is the nominal density
   density.ProjectCoefficient(rho_coeff);

   // * Normalize the density
   density /= rho;

   // * Square the density to get the B matrix
   // density will be needed later, so use a copy
   ParGridFunction bsqrd(density);
   bsqrd *= density; // this is the B matrix; here it is squared
   // coefficient for the mass matrix
   GridFunctionCoefficient bsqrd_coeff(&bsqrd);

   ParBilinearForm *m = new ParBilinearForm(fespace);
   m->AddDomainIntegrator(new MassIntegrator(bsqrd_coeff));
   m->Assemble();
   // shift the eigenvalue corresponding to eliminated dofs to a large value
   m->EliminateEssentialBCDiag(ess_bdr, numeric_limits<double>::min());
   m->Finalize();

   HypreParMatrix *A = a->ParallelAssemble();
   HypreParMatrix *M = m->ParallelAssemble();

#if defined(MFEM_USE_SUPERLU) || defined(MFEM_USE_STRUMPACK)
   Operator * Arow = NULL;
#ifdef MFEM_USE_SUPERLU
   if (slu_solver)
   {
      Arow = new SuperLURowLocMatrix(*A);
   }
#endif
#ifdef MFEM_USE_STRUMPACK
   if (sp_solver)
   {
      Arow = new STRUMPACKRowLocMatrix(*A);
   }
#endif
#endif

   delete a;
   delete m;

   // 8. Define and configure the LOBPCG eigensolver and the BoomerAMG
   //    preconditioner for A to be used within the solver. Set the matrices
   //    which define the generalized eigenproblem A x = lambda M x.
   Solver * precond = NULL;
   if (!slu_solver && !sp_solver && !cpardiso_solver)
   {
      HypreBoomerAMG * amg = new HypreBoomerAMG(*A);
      amg->SetPrintLevel(0);
      precond = amg;
   }
   else
   {
#ifdef MFEM_USE_SUPERLU
      if (slu_solver)
      {
         SuperLUSolver * superlu = new SuperLUSolver(MPI_COMM_WORLD);
         superlu->SetPrintStatistics(false);
         superlu->SetSymmetricPattern(true);
         superlu->SetColumnPermutation(superlu::PARMETIS);
         superlu->SetOperator(*Arow);
         precond = superlu;
      }
#endif
#ifdef MFEM_USE_STRUMPACK
      if (sp_solver)
      {
         STRUMPACKSolver * strumpack = new STRUMPACKSolver(argc, argv, MPI_COMM_WORLD);
         strumpack->SetPrintFactorStatistics(true);
         strumpack->SetPrintSolveStatistics(false);
         strumpack->SetKrylovSolver(strumpack::KrylovSolver::DIRECT);
         strumpack->SetReorderingStrategy(strumpack::ReorderingStrategy::METIS);
         strumpack->DisableMatching();
         strumpack->SetOperator(*Arow);
         strumpack->SetFromCommandLine();
         precond = strumpack;
      }
#endif
#ifdef MFEM_USE_MKL_CPARDISO
      if (cpardiso_solver)
      {
         auto cpardiso = new CPardisoSolver(A->GetComm());
         cpardiso->SetMatrixType(CPardisoSolver::MatType::REAL_STRUCTURE_SYMMETRIC);
         cpardiso->SetPrintLevel(1);
         cpardiso->SetOperator(*A);
         precond = cpardiso;
      }
#endif
   }

   HypreLOBPCG * lobpcg = new HypreLOBPCG(MPI_COMM_WORLD);
   lobpcg->SetNumModes(nev);
   lobpcg->SetRandomSeed(seed);
   lobpcg->SetPreconditioner(*precond);
   lobpcg->SetMaxIter(1000);
   lobpcg->SetTol(1e-8);
   lobpcg->SetPrecondUsageMode(1);
   lobpcg->SetPrintLevel(1);
   lobpcg->SetMassMatrix(*M);
   lobpcg->SetOperator(*A);

   // 9. Compute the eigenmodes and extract the array of eigenvalues. Define a
   //    parallel grid function to represent each of the eigenmodes returned by
   //    the solver.
   Array<double> eigenvalues;
   lobpcg->Solve();
   lobpcg->GetEigenvalues(eigenvalues);
   ParGridFunction x(fespace);

   // TODO : Calculate k_eff using the eigenvalue
   double lambda = eigenvalues[0];
   double k_eff = c2/(lambda + c1);
   double D = sigma_s/(3*sigma_t*sigma_t);
   double k_eff_2ndm = nu_sigma_f/(lambda*D + sigma_a);
   // Print the k_eff
   if (myid == 0)
   {
      cout << "Eigenvalue/ Multiplication factor K_eff = " << k_eff << "." << "\n";
      cout << "K-eff using the 2nd method: " << k_eff_2ndm << "." << "\n";
   }

   /***
    * TODO : Calculate the power distribution
    * * power = (nu_sigma_f/2.4) * phi) * 200 MeV * 1.602e-13 J/MeV * 1e-6 MW/J
    * here phi is the neutron flux and 2.4 is the average nu
    */ 
   ParGridFunction power(fespace);
   ParGridFunction phi(fespace);
   
   // get phi from the first eigenmode
   phi = lobpcg->GetEigenvector(0);

   // get all the grid coefficients
   double all_coeffs_of_power = (nu_sigma_f/2.4) * 200 * 1.602e-13 * 1e-6;

   // project the power coefficient to the power grid function
   power = phi;
   power *= all_coeffs_of_power;
   power *= phi;

   // 10. Save the refined mesh and the modes in parallel. This output can be
   //     viewed later using GLVis: "glvis -np <np> -m mesh -g mode".
   {
      ostringstream mesh_name, mode_name;
      mesh_name << "mesh." << setfill('0') << setw(6) << myid;

      ofstream mesh_ofs(mesh_name.str().c_str());
      mesh_ofs.precision(8);
      pmesh->Print(mesh_ofs);

      // * only need the first eigenmode
      for (int i=0; i<nev; i++)
      {
         // convert eigenvector from HypreParVector to ParGridFunction
         x = lobpcg->GetEigenvector(i);

         mode_name << "mode_" << setfill('0') << setw(2) << i << "."
                   << setfill('0') << setw(6) << myid;

         ofstream mode_ofs(mode_name.str().c_str());
         mode_ofs.precision(8);
         x.Save(mode_ofs);
         mode_name.str("");
      }
   }

   // 11. Send the solution by socket to a GLVis server.
   if (visualization)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      socketstream mode_sock(vishost, visport);
      mode_sock.precision(8);

      for (int i=0; i<nev; i++)
      {
         if ( myid == 0 )
         {
            cout << "Eigenmode " << i+1 << '/' << nev
                 << ", Lambda = " << eigenvalues[i] << endl;
         }

         // convert eigenvector from HypreParVector to ParGridFunction
         x = lobpcg->GetEigenvector(i);

         mode_sock << "parallel " << num_procs << " " << myid << "\n"
                   << "solution\n" << *pmesh << x << flush
                   << "window_title 'Eigenmode " << i+1 << '/' << nev
                   << ", Lambda = " << eigenvalues[i] << "'" << endl;

         char c;
         if (myid == 0)
         {
            cout << "press (q)uit or (c)ontinue --> " << flush;
            cin >> c;
         }
         MPI_Bcast(&c, 1, MPI_CHAR, 0, MPI_COMM_WORLD);

         if (c != 'c')
         {
            break;
         }
      }
      mode_sock.close();
   }

   // 12. Save data in the ParaView format.
   if (paraview)
   {
      //x = lobpcg->GetEigenvector(0);
      ParaViewDataCollection paraview_dc("fluxSol", pmesh);
      paraview_dc.SetPrefixPath("ParaView");
      paraview_dc.SetLevelsOfDetail(order);
      paraview_dc.SetDataFormat(VTKFormat::BINARY);
      paraview_dc.SetHighOrderOutput(true);
      paraview_dc.SetCycle(0);
      paraview_dc.SetTime(0.0);
      paraview_dc.RegisterField("flux", &phi);
      paraview_dc.RegisterField("density", &density);
      paraview_dc.RegisterField("power", &power);
      paraview_dc.Save();
   }



   // 13. Free the used memory.
   delete lobpcg;
   delete precond;
   delete M;
   delete A;
#if defined(MFEM_USE_SUPERLU) || defined(MFEM_USE_STRUMPACK)
   delete Arow;
#endif

   delete fespace;
   if (order > 0)
   {
      delete fec;
   }
   delete pmesh;

   return 0;
}



