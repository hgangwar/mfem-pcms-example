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
using pcms::OmegaHFieldAdapter;
using pcms::MFEMFieldAdapter;

using namespace mfem;

// Function to solve the Thermal problem (Heat Conduction)
// Solves for T_out (Temperature) using U_in (Displacement) from the previous step.
void thermal_solver(
    Mesh &mesh_out
    FiniteElementSpace& T_fespace, 
    const GridFunction& U_in, 
    GridFunction& T_out
) {
    // --- 1. Define Bilinear Form (LHS: Stiffness Matrix) ---
    // A(T, v) = Integral( k * grad(T) . grad(v) ) dOmega
    BilinearForm a_therm(&T_fespace);
    
    // NOTE: In a real coupled problem, 'k' (thermal conductivity) might be a function
    // of the strain derived from U_in, but for simplicity, we use a constant.
    ConstantCoefficient k(50.0); 
    a_therm.AddDomainIntegrator(new DiffusionIntegrator(k));
    a_therm.Assemble();
    
    // --- 2. Define Linear Form (RHS: Heat Source Vector and Boundary Conditions) ---
    // L(v) = Integral( Q * v ) dOmega
    LinearForm b_therm(&T_fespace);
    ConstantCoefficient Q(100.0); // Volumetric heat source
    b_therm.AddDomainIntegrator(new DomainLFIntegrator(Q));
    b_therm.Assemble();
    
    // --- 3. Apply Boundary Conditions and Solve ---
    // Assuming boundary attribute '1' for essential (Dirichlet) BCs
    Array<int> essential_dofs_T;
    T_fespace.Get; 
    
    // Form the linear system: A * X = B
    SparseMatrix A_therm;
    Vector B_therm, X_therm;
    
    // We assume the initial guess T_out is already set (e.g., all zeros or an initial uniform temperature)
    a_therm.FormLinearSystem(essential_dofs_T, T_out, b_therm, A_therm, X_therm, B_therm);
    
    // Solve the system using an iterative solver (e.g., Conjugate Gradient)
    CGSolver cg;
    cg.SetRelTol(1e-6);
    cg.SetMaxIter(500);
    cg.SetPreconditioner(new GSSmoother); 
    cg.Solve(A_therm, B_therm, X_therm);
    
    // Recover the solution (scatter back to the GridFunction T_out)
    a_therm.RecoverFEMSolution(X_therm, b_therm, T_out);
}

// -----------------------------------------------------------------------------

// Function to solve the Mechanical problem (Linear Elasticity)
// Solves for U_out (Displacement) using T_in (Temperature) from the current step.
void mech_solver(
    Mesh& mesh, 
    FiniteElementSpace& U_fespace, 
    const GridFunction& T_in, 
    GridFunction& U_out
) {
    // Material Properties
    double E = 200.0e9;   // Young's Modulus (Pa)
    double nu = 0.3;      // Poisson's Ratio
    double alpha = 1.0e-5; // Coefficient of Thermal Expansion (1/K)
    double T_ref = 293.15; // Reference Temperature (Kelvin)
    
    // Lame parameters for Linear Elasticity
    double lambda = (E * nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
    double mu = E / (2.0 * (1.0 + nu));

    // --- 1. Define Bilinear Form (LHS: Stiffness Matrix) ---
    // A(u, v) = Integral( 2*mu*epsilon(u):epsilon(v) + lambda*div(u)*div(v) ) dOmega
    BilinearForm a_mech(&U_fespace);
    VectorConstantCoefficient lambda_coeff(lambda);
    VectorConstantCoefficient mu_coeff(mu);
    a_mech.AddDomainIntegrator(new ElasticityIntegrator(lambda_coeff, mu_coeff));
    a_mech.Assemble();
    
    // --- 2. Define Linear Form (RHS: Force Vector) ---
    // L(v) = Integral( Body Force . v ) dOmega + Integral( Thermal Force . v ) dOmega
    LinearForm b_mech(&U_fespace);
    
    // a. External Body Force (e.g., Gravity in the last dimension)
    Vector g(mesh.Dimension()); g = 0.0; g(mesh.Dimension() - 1) = -9.81;
    VectorConstantCoefficient gravity(g);
    b_mech.AddDomainIntegrator(new VectorDomainLFIntegrator(gravity));
    
    // b. Thermal Load (The driving coupling term)
    // The force is proportional to the difference between current and reference temperature: (T_in - T_ref)
    
    // 1. Coefficient for (T_in - T_ref)
    ConstantCoefficient T_ref_coeff(T_ref);
    GridFunctionCoefficient T_in_coeff(&T_in);
    SumCoefficient T_diff(1.0, T_in_coeff, -1.0, T_ref_coeff); 
    
    // 2. Coefficient for the factor: (3*lambda + 2*mu) * alpha
    ConstantCoefficient thermo_factor((3.0 * lambda + 2.0 * mu) * alpha);

    // 3. The full coefficient for the thermal load term (product of the two above)
    ProductCoefficient thermal_load_coeff(T_diff, thermo_factor);
    
    // Add the specific thermal load integrator for elasticity
    b_mech.AddDomainIntegrator(new ElasticityThermalLoadIntegrator(thermal_load_coeff));
    b_mech.Assemble();
    
    // --- 3. Apply Boundary Conditions and Solve ---
    // Assuming boundary attribute '2' for essential (Dirichlet) BCs
    Array<int> essential_dofs_U;
    // U_fespace.GetEssentialTrueDofs ...
    
    SparseMatrix A_mech;
    Vector B_mech, X_mech;
    
    // Use the previous solution U_in as the initial guess U_out for the next step
    U_out = U_in; 
    
    a_mech.FormLinearSystem(essential_dofs_U, U_out, b_mech, A_mech, X_mech, B_mech);
    
    // Solve the system using an iterative solver (e.g., MINRES for symmetric systems)
    MINRESSolver minres;
    minres.SetRelTol(1e-6);
    minres.SetMaxIter(500);
    minres.SetPreconditioner(new GSSmoother); 
    minres.Solve(A_mech, B_mech, X_mech);
    
    // Recover the solution
    a_mech.RecoverFEMSolution(X_mech, b_mech, U_out);
}

// -----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
    // 1. Initial Setup (Simplified)
    // For a real problem, you would load a mesh and define the appropriate FE spaces.
    // We'll create a simple 2D mesh for demonstration.
    Mesh mesh(10, 10, Element::QUADRILATERAL);
    int order = 1;

    // 2. Define Finite Element Spaces
    H1_FECollection T_fec(order, mesh.Dimension());
    FiniteElementSpace T_fespace(&mesh, &T_fec);
    
    H1_FECollection U_fec(order, mesh.Dimension());
    FiniteElementSpace U_fespace(&mesh, &U_fec, mesh.Dimension()); // Vector space for displacement

    std::cout << "Number of T unknowns: " << T_fespace.Get=VSize() << std::endl;
    std::cout << "Number of U unknowns: " << U_fespace.Get=VSize() << std::endl;
    
    // 3. Define Solution/Storage Vectors
    // T_k_minus_1 and U_k_minus_1 hold the solution from the *previous* iteration.
    // T_k and U_k hold the solution from the *current* iteration.
    GridFunction T_k_minus_1(&T_fespace), T_k(&T_fespace);
    GridFunction U_k_minus_1(&U_fespace), U_k(&U_fespace);

    // 4. Initial Guess and Convergence Parameters
    T_k_minus_1 = 300.0; // Start with a uniform temperature of 300 K
    U_k_minus_1 = 0.0;   // Start with zero displacement
    
    int max_iter = 50;
    double tolerance = 1e-6;
    double T_diff = tolerance * 2.0; // Initialize high to enter the loop
    double U_diff = tolerance * 2.0; 

    // 5. Staggered Iteration Loop (Gauss-Seidel-like)
    std::cout << "\nStarting Staggered Thermo-Mechanical Iteration..." << std::endl;

    for (int k = 1; k <= max_iter; ++k) {
        
        // --- STEP 1: SOLVE THERMAL ---
        // T_k = Solver(U_{k-1}) - Thermal solution depends on mechanical strain
        thermal_solver(T_fespace, U_k_minus_1, T_k);
        
        // --- STEP 2: SOLVE MECHANICAL ---
        // U_k = Solver(T_k) - Mechanical displacement depends on current temperature
        mech_solver(mesh, U_fespace, T_k, U_k);

        // --- STEP 3: CHECK CONVERGENCE ---
        T_diff = T_k.DistanceTo(T_k_minus_1);
        U_diff = U_k.DistanceTo(U_k_minus_1);

        std::cout << "Iter " << k 
                  << ": T_res = " << T_diff 
                  << ", U_res = " << U_diff << std::endl;
        
        // Update previous solutions for the next iteration
        T_k_minus_1 = T_k;
        U_k_minus_1 = U_k;

        if (T_diff < tolerance && U_diff < tolerance) {
            std::cout << "\nCONVERGED in " << k << " iterations. 🎉" << std::endl;
            break;
        }
        
        if (k == max_iter) {
            std::cout << "\nWARNING: Max iterations reached. Solution did not converge." << std::endl;
        }
    }

    // 6. Visualization (Optional)
    // char vishost[] = "localhost";
    // int visport = 19916;
    // socketstream sout_T(vishost, visport);
    // sout_T << "solution\n" << mesh << T_k << "caption 'Final Temperature T'" << std::endl;
    // socketstream sout_U(vishost, visport);
    // sout_U << "solution\n" << mesh << U_k << "caption 'Final Displacement U'" << std::endl;

    return 0;
}