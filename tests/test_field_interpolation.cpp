
#include <Omega_h_mesh.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_for.hpp>
#include <pcms/transfer_field2.h>
#include "pcms/omega_h_field.h"
#include "pcms/omega_h_field2.h"
#include "pcms/create_field.h"
#include <Kokkos_Core.hpp>
#include <vector>

using pcms::Real;

int main()
{
  auto lib = Omega_h::Library{};
  auto world = lib.world();
  auto mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 0, 100, 100, 0, false);
  auto layout =
    pcms::CreateLagrangeLayout(mesh, 1, 1, pcms::CoordinateSystem::Cartesian);
  const auto nverts = mesh.nents(0);
  auto mesh_coords = mesh.coords();
  auto f = [](Real x, Real y) { return -0.3 * x + 0.5 * y; };
  Omega_h::Write<Real> test_f(nverts);
  Omega_h::parallel_for(
    nverts, OMEGA_H_LAMBDA(int i) {
      Real x = mesh_coords[2 * i + 0];
      Real y = mesh_coords[2 * i + 1];
      test_f[i] = f(x, y);
    });
  auto field = layout->CreateField();
  auto interpolated = layout->CreateField();
  field->SetDOFHolderData(pcms::make_const_array_view(test_f));

  pcms::interpolate_field2(*field, *interpolated);
  auto interpolated_dof = interpolated->GetDOFHolderData();
  auto original_dof = field->GetDOFHolderData();
  assert(interpolated_dof.size() == original_dof.size());

  // assumes that GetDOFHolderData will return a host view
  for (std::size_t i = 0; i < interpolated_dof.size(); ++i) {
    const double a = interpolated_dof[i];
    const double b = original_dof[i];

    const bool within_rel =
        std::abs(a - b) <= 0.001 * std::abs(b);

    const bool within_abs =
        std::abs(a - b) <= 1e-10;

    assert(within_rel || within_abs);
  }
  return 0;
}

