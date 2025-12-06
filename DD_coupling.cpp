#include "mfem.hpp"
#include <iostream>
#include <fstream>

using namespace mfem;
using namespace std;

// ------------------------------------------------------------
// Assign element attributes based on x-coordinate:
//   attribute = 1 if x ≤ 0.6  (active region)
//   attribute = 2 otherwise   (inactive region)
// ------------------------------------------------------------
void AssignAttributesByX(Mesh &mesh)
{
    for (int e = 0; e < mesh.GetNE(); e++)
    {
        Vector center;
        mesh.GetElementCenter(e, center);
        double x = center[0];

        if (x <= 0.6)
            mesh.SetAttribute(e, 1);
        else
            mesh.SetAttribute(e, 2);
    }
}

// ------------------------------------------------------------
// Select DOFs that lie on the interior plane x = 0.5
// Produces an essential-dof list: ess_dofs
// ------------------------------------------------------------
void MarkInteriorPlaneDOFs(const Mesh &mesh,
                           const FiniteElementSpace &fes,
                           Array<int> &ess_dofs)
{
  const double xc = 0.5;
  const double tol = 1e-12;

  Array<int> vdof_marker(fes.GetVSize());
  vdof_marker = 0;

  const int dim = mesh.Dimension();

  for (int v = 0; v < mesh.GetNV(); v++)
  {
    const double *vx = mesh.GetVertex(v);  // <-- const-safe method

    double x = vx[0];                      // X coordinate

    if (fabs(x - xc) < tol)
    {
      Array<int> vdofs;
      fes.GetVertexDofs(v, vdofs);
      for (int j = 0; j < vdofs.Size(); j++)
        vdof_marker[vdofs[j]] = 1;
    }
  }

  fes.GetEssentialVDofs(vdof_marker, ess_dofs);
}


int main(int argc, char *argv[])
{
    std::string mesh_file = "$HOME/src/mfem-pcms-example/mesh/cube.msh";

    Mesh mesh = new mfem::Mesh(mesh_file.c_str(), 1, 1, true);

    // Assign attributes region-wise
    AssignAttributesByX(mesh);

    // FE space
    int order = 1;
    H1_FECollection fec(order, mesh.Dimension());
    FiniteElementSpace fes(&mesh, &fec);

    GridFunction sol(&fes);
    sol = 0.0;

    ConstantCoefficient one(1.0);

    // ------------------------------------------------------------
    // Attribute mask (region 1 active)
    // ------------------------------------------------------------
    int max_attr = mesh.attributes.Max();
    Array<int> attr_mask(max_attr);
    attr_mask = 0;
    attr_mask[0] = 1;
    attr_mask[1] = 1;// attribute 1 active

    // ------------------------------------------------------------
    // Bilinear form (masked)
    // ------------------------------------------------------------
    BilinearForm a(&fes);
    a.SetAssemblyLevel(AssemblyLevel::FULL);
    a.AddDomainIntegrator(new DiffusionIntegrator(one), attr_mask);
    a.Assemble();

    // RHS (masked)
    LinearForm b(&fes);
    b.AddDomainIntegrator(new DomainLFIntegrator(one), attr_mask);
    b.Assemble();

    // ------------------------------------------------------------
    // INTERNAL essential BC at x = 0.5
    // ------------------------------------------------------------
    Array<int> interior_ess_dofs;
    MarkInteriorPlaneDOFs(mesh, fes, interior_ess_dofs);

    // Inhomogeneous BC value (0.0)
    sol = 0.0;

    // ------------------------------------------------------------
    // Solve masked system with internal BC
    // ------------------------------------------------------------
    SparseMatrix A;
    Vector B, X;
    a.FormLinearSystem(interior_ess_dofs, sol, b, A, X, B);

    GSSmoother M(A);
    PCG(A, M, B, X, 1, 300, 1e-12, 0.0);

    a.RecoverFEMSolution(X, b, sol);

    // ------------------------------------------------------------
    // Output solution
    // ------------------------------------------------------------
    {
        ofstream vtk("solution.vtk");
        sol.SaveVTK(vtk, "solution", 1);
    }

    // ------------------------------------------------------------
    // Output mask field (region 1 indicator)
    // ------------------------------------------------------------
    GridFunction mask_gf(&fes);
    mask_gf = 0.0;

    for (int e = 0; e < mesh.GetNE(); e++)
    {
        int attr = mesh.GetAttribute(e);
        if (attr == 1)
        {
            Array<int> vdofs;
            fes.GetElementVDofs(e, vdofs);
            for (int j = 0; j < vdofs.Size(); j++)
                mask_gf(vdofs[j]) = 1.0;
        }
    }

    {
        ofstream mvtk("mask.vtk");
        mask_gf.SaveVTK(mvtk, "mask", 1);
    }

    cout << "Wrote: solution.vtk and mask.vtk\n";
    return 0;
}
