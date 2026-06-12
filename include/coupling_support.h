#pragma once

#include "mfem.hpp"
#include <pcms/pcms.h>
#include <pcms/create_field.h>
#include <Omega_h_array.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_mesh.hpp>
#include "pcms/coupler2.h"
namespace support
{

using dtype = pcms::Real;

// -----------------------------
// FEMSystem
// -----------------------------
struct FEMSystem
{
  mfem::ParMesh* pmesh = nullptr;
  mfem::H1_FECollection* fec = nullptr;
  mfem::ParFiniteElementSpace* fes = nullptr;

  mfem::ParBilinearForm* a = nullptr;
  mfem::ParLinearForm* b = nullptr;
  mfem::ParGridFunction* x = nullptr;

  mfem::Array<int> ess_tdofs; // essential TRUE dofs for elimination
};

// -----------------------------
// System Parameters
// -----------------------------
struct ThermalParams
{
  std::array<double, 3> size; // lx, ly, lz
  std::array<int, 3> ne;      // nx, ny, nz
  double q_total;             // total heat [W]
  double kappa;               // thermal conductivity
  double rho;                 // density
  double cp;                  // heat capacity (unused here, steady)
  double h_flux;              // peak heat transfer coefficient for flux BC
  double h_conv;              // convection coefficient
  double T_conv;              // ambient temperature for convection BC
  double T_dirichlet;         // Dirichlet boundary temperature
};

// -----------------------------
// Coupling metadata
// -----------------------------
struct Coupling
{
  std::string name;                                  // Coupler name
  std::unique_ptr<pcms::Coupler2> cpl;                // Coupler instance
  std::vector<std::string> app_names;                // Attached apps
  std::vector<std::string> field_names;              // Attached fields
  std::map<std::string, pcms::Application2*> apps;    // App name -> pointer
  std::map<std::string, pcms::CoupledField*> fields; // App field -> pointer
  bool isServer = false;
};

// -----------------------------
// Exact solution coefficient
// -----------------------------
class ExactTempCoeff : public mfem::Coefficient
{
public:
  double Eval(mfem::ElementTransformation& T,
              const mfem::IntegrationPoint& ip) override;
};

// -----------------------------
// ParaView output helpers
// -----------------------------
struct OutputPack
{
  mfem::ParaViewDataCollection pvd;
  mfem::ParGridFunction exact;
  mfem::ParGridFunction err;

  OutputPack(const std::string& collection, mfem::ParMesh& pm,
             mfem::ParFiniteElementSpace& fes);
};

// -----------------------------
// FEM functions
// -----------------------------
FEMSystem Init_FEMSystem(mfem::ParMesh* pmesh, int order, double kappa_val);

void DestroyFEMSystem(FEMSystem& sys);

long SolveSystem(FEMSystem& sys, const std::string& solver_type,
                 const std::string& prec_type, double rel_tol, int max_iter);

// -----------------------------
// Utilities
// -----------------------------
double DefaultTolX(const mfem::Mesh& mesh);

void shift_meshX(Omega_h::Mesh& mesh, double dx);

Omega_h::Write<Omega_h::I8> create_mask(Omega_h::Mesh& mesh,
                                        const char* tag_name, int tag_value);

double RMSDiff(const std::vector<std::pair<double, double>>& a,
               const std::vector<std::pair<double, double>>& b);

long ComputeRMS(const Omega_h::Read<dtype>& a, const Omega_h::Read<dtype>& b);

std::vector<std::pair<double, double>> ExtractVertexLineTrace(
  const Omega_h::Mesh& mesh, Omega_h::Read<Omega_h::Real> field_v, double xline,
  double tol);

void FillTagOnXLineFromTrace(
  Omega_h::Mesh& mesh, const std::vector<std::pair<double, double>>& trace,
  double x_line, double tol, const char* tag_name = "temp");

void ApplyBoundaryTraceByAttr(
  mfem::ParMesh& pmesh, mfem::ParGridFunction& gf, int bdr_attr,
  const std::vector<std::pair<double, double>>& trace, double tol);

void ApplyBoundaryConstantByAttr(mfem::ParMesh& pmesh,
                                 mfem::ParGridFunction& gf, int bdr_attr,
                                 double value);

void ReportBdrAttrStats(const mfem::ParMesh& pmesh,
                        const mfem::ParGridFunction& T, int bdr_attr,
                        const char* name);

void ReportTraceStats(const std::vector<std::pair<double, double>>& tr,
                      const char* name);

std::vector<std::pair<double, double>> RelaxTrace(
  const std::vector<std::pair<double, double>>& old_t,
  const std::vector<std::pair<double, double>>& new_t, double omega);

void SaveFields(OutputPack& out, const FEMSystem& sys, int it);

void write_oh_mesh(Omega_h::Mesh& mesh, const std::string& path);

mfem::Mesh read_mfem_mesh(const std::string& path);

void reset_mfem_attributes(mfem::Mesh& mesh, int attr = 1);

void remove_oh_tag(Omega_h::Mesh& mesh, const std::string& name);

void SaveParaview(mfem::ParMesh& pmesh, mfem::ParGridFunction& x,
                  const std::string& collection = "thermal_solution",
                  const std::string& field_name = "Temperature", int cycle = 0,
                  double time = 0.0);

//--------------------------------------------------------------
// Init_Coupler<TAdapter>
//
// Templates must stay in the header unless explicitly instantiated
// in the .cpp for every adapter type you use.
//--------------------------------------------------------------
template <typename Adapter_type>
Coupling Init_Coupler(MPI_Comm comm, const std::string& name,
                      const std::vector<std::string>& app_names,
                      const std::vector<std::string>& field_names,
                      bool isServer, const redev::Partition ptn,
                      Adapter_type& Adapter)
{
  Coupling cp;
  cp.name = name;
  cp.app_names = app_names;
  cp.field_names = field_names;
  cp.isServer = isServer;

  if (app_names.size() != field_names.size()) {
    throw std::runtime_error(
      "Mismatch: app_names and field_names must be of the same size.");
  }

  cp.cpl = std::make_unique<pcms::Coupler2>(name, comm, isServer, ptn);

  for (size_t i = 0; i < app_names.size(); ++i) {
    auto* app = cp.cpl->AddApplication(app_names[i]);


    cp.fields[app_names[i]] = app->AddField(field_names[i], std::move(Adapter));
    cp.apps[app_names[i]] = app;
  }

  return cp;
}
void initializeFieldWithGids(pcms::FieldT<pcms::Real>* field,
                             pcms::Real multiplier);
template <typename dtype>
Coupling Init_Coupler_OH(MPI_Comm comm, const std::string& name,
                         const std::vector<std::string>& app_names,
                         const std::vector<std::string>& field_names,
                         bool isServer, const redev::Partition ptn,
                         Omega_h::Mesh& mesh,
                         const Omega_h::Write<Omega_h::I8> is_overlap)
{
  if (app_names.size() != field_names.size())
    throw std::runtime_error("app_names and field_names size mismatch.");

  Coupling cp;
  cp.name = name;
  cp.app_names = app_names;
  cp.field_names = field_names;
  cp.isServer = isServer;

  cp.cpl = std::make_unique<pcms::Coupler2>(name, comm, isServer, ptn);

  for (size_t i = 0; i < app_names.size(); ++i) {
    auto* app = cp.cpl->AddApplication(app_names[i]);
    auto& layout = app->AddLayout(
        "temp",
        pcms::CreateLagrangeLayout(mesh, 1, 1, pcms::CoordinateSystem::Cartesian));
    auto field = layout.CreateFieldReal();
    auto* field_ptr = field.get();
    initializeFieldWithGids(field_ptr, 0.0);
    app->AddField(field_names[i], std::move(field));

    cp.apps[app_names[i]] = app;
  }

  return cp;
}
double ComputeAbsoluteError(const Omega_h::Mesh& mesh);
void PrintTempStats(const Omega_h::Mesh& mesh, const std::string& name,
                    int itr);

double ComputeAbsoluteError(const mfem::ParMesh& pmesh,
                            const mfem::ParGridFunction& x);

void PrintTempStats(const mfem::ParMesh& pmesh, const mfem::ParGridFunction& x,
                    const std::string& name, int itr);

} // namespace support