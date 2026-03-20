#include "mfem.hpp"
#include <iostream>
#include <memory>
#include <cmath>

using namespace mfem;

// Exact scalar field to transfer.
// Linear so it is represented exactly with H1 order 1 on this mesh.
double exact_function(const Vector &x)
{
   return x(0) + 2.0 * x(1);
}

int main(int argc, char *argv[])
{
   // ------------------------------------------------------------
   // 1) Build a simple parent mesh: [0,1] x [0,1], 4 x 2 quads
   // ------------------------------------------------------------
   Mesh mesh = Mesh::MakeCartesian2D(4, 2, Element::QUADRILATERAL,
                                     true, 1.0, 1.0);

   // Mark left half as attr=1, right half as attr=2
   for (int e = 0; e < mesh.GetNE(); ++e)
   {
      Array<int> verts;
      mesh.GetElementVertices(e, verts);

      double xc = 0.0;
      for (int i = 0; i < verts.Size(); ++i)
      {
         const double *v = mesh.GetVertex(verts[i]);
         xc += v[0];
      }
      xc /= verts.Size();

      mesh.GetElement(e)->SetAttribute((xc < 0.5) ? 1 : 2);
   }

   std::cout << "Parent mesh elements: " << mesh.GetNE() << "\n";
   std::cout << "Parent mesh attributes:";
   for (int i = 0; i < mesh.attributes.Size(); ++i)
   {
      std::cout << " " << mesh.attributes[i];
   }
   std::cout << "\n";

   // ------------------------------------------------------------
   // 2) Create a submesh from domain attribute {2}
   // ------------------------------------------------------------
   Array<int> attrs(1);
   attrs[0] = 2;

   SubMesh submesh = SubMesh::CreateFromDomain(mesh, attrs);

   std::cout << "Submesh elements: " << submesh.GetNE() << "\n";

   // Print parent->submesh element map for sanity
   const Array<int> &parent_elem_ids = submesh.GetParentElementIDMap();
   std::cout << "Submesh element -> parent element map:\n";
   for (int se = 0; se < parent_elem_ids.Size(); ++se)
   {
      std::cout << "  sub elem " << se
                << " -> parent elem " << parent_elem_ids[se] << "\n";
   }

   // ------------------------------------------------------------
   // 3) Build FE spaces
   // ------------------------------------------------------------
   const int order = 1;
   H1_FECollection fec(order, mesh.Dimension());

   FiniteElementSpace parent_fes(&mesh, &fec);
   FiniteElementSpace sub_fes(&submesh, &fec);

   GridFunction parent_gf(&parent_fes);
   GridFunction sub_gf(&sub_fes);
   GridFunction sub_exact(&sub_fes);

   // ------------------------------------------------------------
   // 4) Project exact field on parent and submesh separately
   // ------------------------------------------------------------
   FunctionCoefficient exact_coeff(exact_function);

   parent_gf.ProjectCoefficient(exact_coeff);
   sub_exact.ProjectCoefficient(exact_coeff);

   // Initialize destination to a recognizable junk value
   sub_gf = -999.0;

   // ------------------------------------------------------------
   // 5) Create and apply transfer map: parent -> submesh
   // ------------------------------------------------------------
   auto transfer_map = SubMesh::CreateTransferMap(parent_gf, sub_gf);
   transfer_map.Transfer(parent_gf, sub_gf);

   // MFEM also has SubMesh::Transfer(parent_gf, sub_gf), but here
   // we explicitly create the map once and reuse it.

   // ------------------------------------------------------------
   // 6) Verify correctness
   // ------------------------------------------------------------
   GridFunction diff(&sub_fes);
   diff = sub_gf;
   diff -= sub_exact;

   const double l2_err   = diff.Norml2();
   const double linf_err = diff.Normlinf();

   std::cout << "Transfer verification on submesh:\n";
   std::cout << "  ||transferred - exact||_L2   = " << l2_err << "\n";
   std::cout << "  ||transferred - exact||_Linf = " << linf_err << "\n";

   const double tol = 1e-12;
   if (linf_err < tol)
   {
      std::cout << "PASS: transfer map is correct.\n";
      return 0;
   }
   else
   {
      std::cout << "FAIL: transfer map is incorrect.\n";
      return 1;
   }
}
