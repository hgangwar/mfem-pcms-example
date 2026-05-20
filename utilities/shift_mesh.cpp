#include <mfem.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>

#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include "schwartz_coupling_support.h"
using namespace mfem;

// --------------------------------------------------
// Helper: strip extension
std::string strip_ext(const std::string &name)
{
    size_t pos = name.find_last_of('.');
    if (pos == std::string::npos) { return name; }
    return name.substr(0, pos);
}



void write_mfem_mesh(const mfem::Mesh &mesh, const std::string &path)
{
    std::ofstream os(path);
    if (!os)
    {
        throw std::runtime_error("Failed to open MFEM output file: " + path);
    }
    mesh.Print(os);
}

// --------------------------------------------------
// Omega_h read/write
Omega_h::Mesh read_oh_mesh(Omega_h::Library &lib, const std::string &path)
{
    Omega_h::Mesh mesh(&lib);
    Omega_h::binary::read(path, lib.world(), &mesh);
    return mesh;
}



// --------------------------------------------------
int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc < 3)
    {
        if (rank == 0)
        {
            std::cerr << "Usage:\n"
                      << argv[0] << " input.mesh input.osh [shift]\n";
        }
        MPI_Finalize();
        return 1;
    }

    const std::string mfem_in = argv[1];
    const std::string oh_in   = argv[2];

    double shift = 0.4;
    if (argc > 3)
    {
        shift = std::atof(argv[3]);
    }

    const std::string mfem_out = strip_ext(mfem_in) + "_shifted.mesh";
    const std::string oh_out   = strip_ext(oh_in)   + "_shifted.osh";

    {
        Omega_h::Library lib(&argc, &argv);

        // ============================================================
        // MFEM
        // ============================================================
        //mfem::Mesh mfem_mesh = support::read_mfem_mesh(mfem_in);
      mfem::Mesh mfem_mesh(mfem_in.c_str(), 1, 1);
        const int nverts_mfem = mfem_mesh.GetNV();
        for (int i = 0; i < nverts_mfem; ++i)
        {
            double *v = mfem_mesh.GetVertex(i);
            v[0] += shift;
        }

        write_mfem_mesh(mfem_mesh, mfem_out);

        if (rank == 0)
        {
            std::cout << "[MFEM] Wrote: " << mfem_out << "\n";
        }

        // ============================================================
        // Omega_h
        // ============================================================
        Omega_h::Mesh oh_mesh = read_oh_mesh(lib, oh_in);

        const int dim    = oh_mesh.dim();
        const int nverts = oh_mesh.nverts();

        auto coords = oh_mesh.coords();
        Omega_h::Write<Omega_h::Real> new_coords(coords.size());

        for (int i = 0; i < nverts; ++i)
        {
            new_coords[dim * i + 0] = coords[dim * i + 0] + shift;
            for (int d = 1; d < dim; ++d)
            {
                new_coords[dim * i + d] = coords[dim * i + d];
            }
        }

        Omega_h::Mesh shifted = oh_mesh;
        shifted.set_coords(Omega_h::Reals(new_coords));

        support::write_oh_mesh(shifted, oh_out);

        if (rank == 0)
        {
            std::cout << "[Omega_h] Wrote: " << oh_out << "\n";
        }
    }

    MPI_Finalize();
    return 0;
}