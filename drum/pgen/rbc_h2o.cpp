/** @file rayleigh_benard.cpp
 * @brief Rayleigh-Benard convection in planetary atmospheres
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Tuesday May 22, 2021 08:07:22 PDT
 * @bug No known bugs.
 */

// C/C++ header
#include <ctime>
#include <sstream>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../math/interpolation.h"
#include "../globals.hpp"
#include "../utils/utils.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "../physics/physics.hpp"

namespace math {
  #include "../math/core.h"
};

// molecules
enum {ivapor = 1};

// global parameters
Real grav, P0, T0, Tmin, radius, omega;
bool use_polar_beta;

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  AllocateUserOutputVariables(4+NVAPOR);
  SetUserOutputVariableName(0, "temp", "temperature", "K");
  SetUserOutputVariableName(1, "theta", "potential temperature", "K");
  SetUserOutputVariableName(2, "thetav", "virtual potential temperature", "K");
  SetUserOutputVariableName(3, "mse", "moist static energy", "J/kg");
  for (int n = 0; n < NVAPOR; ++n)
    SetUserOutputVariableName(4+i, "rh" + std::tostring(n+1), "relative humidity");
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
{
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        user_out_var(0,k,j,i) = pthermo->GetTemp(phydro->w.at(k,j,i));
        user_out_var(1,k,j,i) = PotentialTemp(phydro->w.at(k,j,i), P0, pthermo);
        // theta_v
        user_out_var(2,k,j,i) = user_out_var(1,k,j,i)*pthermo->RovRd(phydro->w.at(k,j,i));
        // mse
        user_out_var(3,k,j,i) = MoistStaticEnergy(phydro->w, grav*pcoord->x1v(i),
          pthermo, ppart, k, j, i);
        for (int n = 0; n < NVAPOR; ++n)
          user_out_var(4+n,k,j,i) = RelativeHumidity(phydro->w.at(k,j,i), ivapor+n, pthermo);
      }
}

void Forcing(MeshBlock *pmb, Real const time, Real const dt,
    AthenaArray<Real> const &w, AthenaArray<Real> const &bcc, AthenaArray<Real> &u)
{
  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        if (use_polar_beta) {
          Real x2 = pmb->pcoord->x2v(j);
          Real x3 = pmb->pcoord->x3v(k);
          Real dist = sqrt(x2*x2 + x3*x3);

          Real fcor = -omega*math::sqr(dist/radius);
          u(IM2,k,j,i) += dt*fcor*w(IDN,k,j,i)*w(IM3,k,j,i);
          u(IM3,k,j,i) -= dt*fcor*w(IDN,k,j,i)*w(IM2,k,j,i);
        }
      }
}

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  grav = - pin->GetReal("hydro", "grav_acc1");
  P0 = pin->GetReal("problem", "P0");
  T0 = pin->GetReal("problem", "T0");
  use_polar_beta = pin->GetOrAddBoolean("problem", "use_polar_beta", false);
  if (use_polar_beta) {
    radius = pin->GetReal("problem", "radius");
    omega = pin->GetReal("hydro", "OmegaZ");
  }
  Tmin = pin->GetReal("hydro", "min_tem");

  EnrollUserExplicitSourceFunction(Forcing);
}

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  std::stringstream msg;
  Real gamma = pin->GetReal("hydro", "gamma");

  // construct a 1D pseudo-moist adiabat with given relative humidity
  Real x1min = pmy_mesh->mesh_size.x1min;
  Real x1max = pmy_mesh->mesh_size.x1max;
  Real x1rat = pmy_mesh->mesh_size.x1rat;

  Real dz, **w1, *z1, *p1, *t1;
  if (x1rat != 1.0) {
    dz = (x1max - x1min)*(x1rat - 1.)/(pow(x1rat, pmy_mesh->mesh_size.nx1) - 1.);
    dz = std::min(dz, dz*pow(x1rat, pmy_mesh->mesh_size.nx1))/2.;
  } else {
    dz = (x1max - x1min)/pmy_mesh->mesh_size.nx1/2.;
  }
  int nx1 = (int)((x1max - x1min)/dz);
  NewCArray(w1, nx1, NHYDRO+2);
  z1 = new Real [nx1];
  p1 = new Real [nx1];
  t1 = new Real [nx1];

  // 1.1 estimate surface temperature and pressure
  Real Rd = pthermo->GetRd();
  Real cp = gamma/(gamma - 1.)*Rd;
  Real Ts = T0 - grav/cp*x1min;
  Real Ps = P0*pow(Ts/T0, cp/Rd);
  int max_iter = 200, iter = 0;

  for (int n = 0; n < NVAPOR; ++n) {
    Real qv = pin->GetReal("problem", "qvapor" + std::tostring(n+1))/1.E3;
    w1[0][ivapor+n] = qv;
  }
  z1[0] = x1min;
  for (int i = 1; i < nx1; ++i)
    z1[i] = z1[i-1] + dz;

  Real t0, p0;
  while (iter++ < max_iter) {
    pthermo->ConstructAtmosphere(w1, Ts, Ps, grav, dz, nx1, Adiabat::pseudo);

    // 1.2 replace adiabatic atmosphere with isothermal atmosphere if temperature is too low
    int ii = 0;
    for (; ii < nx1; ++ii)
      if (pthermo->Temp(w1[ii]) < Tmin) break;
    Real Tv = w1[ii][IPR]/(w1[ii][IDN]*Rd);
    for (int i = ii; i < nx1; ++i) {
      w1[i][IPR] = w1[ii][IPR]*exp(-grav*(z1[i] - z1[ii])/(Rd*Tv));
      w1[i][IDN] = w1[i][IPR]/(Rd*Tv);
      for (int n = 1; n < NMASS; ++n)
        w1[i][n] = w1[ii][n];
    }

    // 1.3 find TP at z = 0
    for (int i = 0; i < nx1; ++i) {
      p1[i] = w1[i][IPR];
      t1[i] = pthermo->Temp(w1[i]);
    }
    p0 = interp1(0., p1, z1, nx1);
    t0 = interp1(0., t1, z1, nx1);

    Ts += T0 - t0;
    Ps *= P0/p0;
    if ((fabs(T0 - t0) < 0.01) && (fabs(P0/p0 - 1.) < 1.E-4)) break;
  }

  if (iter > max_iter) {
    msg << "### FATAL ERROR in problem generator"
        << std::endl << "maximum iteration reached."
        << std::endl << "T0 = " << t0
        << std::endl << "P0 = " << p0;
    ATHENA_ERROR(msg);
  }

  // setup initial condition
  //srand(Globals::my_rank + time(0));
  for (int i = is; i <= ie; ++i) {
    Real buf[NHYDRO+2*NVAPOR];
    interpn(buf, &pcoord->x1v(i), *w1, z1, &nx1, 1, NHYDRO+2*NVAPOR);
    buf[IVX] = buf[IVY] = buf[IVZ] = 0.;

    // set gas concentration
    for (int n = 0; n < NHYDRO; ++n)
      for (int k = ks; k <= ke; ++k)
        for (int j = js; j <= je; ++j)
          phydro->w(n,k,j,i) = buf[n];

    // set cloud concentration
    Particle *pp = ppart;
    for (int n = 0; n < NVAPOR; ++n) {
      if (pp == nullptr) {
        msg << "### FATAL ERROR in problem generator"
            << std::endl << "Vapor #" << n
            << "does not associate with any particle.";
        ATHENA_ERROR(msg);
      }
      for (int k = ks; k <= ke; ++k)
        for (int j = js; j <= je; ++j)
          pp->c(0,k,j,i) = buf[NHYDRO+n] + buf[NHYDRO+NVAPOR+n];
      pp = pp->next;
    }
  }

  ppart->Initialize();
  pphy->Initialize(phydro->w);
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);

  FreeCArray(w1);
  delete[] z1;
  delete[] p1;
  delete[] t1;
}
