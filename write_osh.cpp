//
// Created by gangwh on 2/9/26.
//
// write_box_osh.cpp
#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_file.hpp>   // binary::write

#include <mpi.h>

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  {
    int nx = 30, ny = 30;

    Omega_h::Library lib(&argc, &argv);

    // Same geometry as MFEM MakeCartesian2D(nx,ny,TRIANGLE,true,0.6,1.0)
    Omega_h::Mesh mesh = Omega_h::build_box(
      lib.world(),
      OMEGA_H_SIMPLEX,
      0.6, 1.0, 0.0,
      nx,  ny,  0
    );

    // Write .osh (single file in serial; in MPI it’s still valid to call on all ranks)
    Omega_h::binary::write("box_tri.osh", &mesh);
  }
  MPI_Finalize();
  return 0;
}
