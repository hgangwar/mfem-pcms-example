//
// Created by gangwh on 2/9/26.
//
// write_box_osh.cpp
#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_file.hpp>
#include <mpi.h>

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  {
    int nx = 30, ny = 30;

    Omega_h::Library lib(&argc, &argv);

    Omega_h::Mesh mesh = Omega_h::build_box(
      lib.world(),
      OMEGA_H_SIMPLEX,
      0.1, 1.0, 0.0,
      nx,  ny,  0
    );

    // Write .osh
    Omega_h::binary::write("cube.osh", &mesh);
  }
  MPI_Finalize();
  return 0;
}
