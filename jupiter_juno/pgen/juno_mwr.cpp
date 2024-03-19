/* -------------------------------------------------------------------------------------
 * Juno MWR inversion Program
 *
 * Contributer:
 * Cheng Li, University of Michigan
 *
 * Year: 2022
 * Contact: chengcli@umich.edu
 * Reference: TBD
 * -------------------------------------------------------------------------------------
 */

// @sect3{Include files}

// C/C++ header
#include <iostream>
#include <algorithm>

// Athena++ headers
// These input files are the same as those in the @ref straka. No additional comments
// are needed.
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../thermodynamics/thermodynamic_funcs.hpp"

// Here we include more input files
#include "../thermodynamics/molecules.hpp"
#include "../thermodynamics/thermodynamic_funcs.hpp"
#include "../thermodynamics/vapors/sodium_vapors.hpp"
#include "../globals.hpp"
#include "../utils/utils.hpp"
#include "../math/interpolation.h"
#include "../math/linalg.h"
#include "../debugger/debugger.hpp"
#include "../radiation/radiation.hpp"
#include "../radiation/microwave/mwr_absorbers.hpp"
#include "../inversion/profile_inversion.hpp"
#include "../planet_data/jup_fletcher16_cirs.hpp"
#include "../planet_data/planets.hpp" // centric2graphic
#include "../scalars/scalars.hpp"

// @sect3{Preamble}

// molecules
enum {iH2O = 1, iNH3 = 2};
enum {ion = 0, iNa = 1, iKCl = 2};
Real grav, P0, T0, Tmin, xHe, xCH4, clat;
Real xH2S, xNa, xKCl, metallicity;

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  AllocateUserOutputVariables(4+NVAPOR);
  SetUserOutputVariableName(0, "temp", "temperature", "K");
  SetUserOutputVariableName(1, "theta", "potential temperature", "K");
  SetUserOutputVariableName(2, "thetav", "virtual potential temperature", "K");
  SetUserOutputVariableName(3, "mse", "moist static energy", "J/kg");
  for (int n = 1; n <= NVAPOR; ++n) {
    std::string name = "rh" + std::to_string(n);
    SetUserOutputVariableName(3+n, name.c_str(), "relative humidity");
  }

  // define inversion
  pfit = new ProfileInversion(this, pin);
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
{
  Real dq[1+NVAPOR], rh;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        user_out_var(0,k,j,i) = pthermo->GetTemp(phydro->w.at(k,j,i));
        user_out_var(1,k,j,i) = PotentialTemp(phydro->w.at(k,j,i), P0, pthermo);
        // theta_v
        user_out_var(2,k,j,i) = user_out_var(1,k,j,i)*pthermo->RovRd(phydro->w.at(k,j,i));
        // mse
        user_out_var(3,k,j,i) = MoistStaticEnergy(phydro->w.at(k,j,i), grav*pcoord->x1v(i), pthermo);
        for (int n = 1; n <= NVAPOR; ++n)
          user_out_var(3+n,k,j,i) = RelativeHumidity(phydro->w.at(k,j,i), n, pthermo);
      }
}

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  grav = - pin->GetReal("hydro", "grav_acc1");
  P0 = pin->GetReal("hydro", "reference_pressure");
  T0 = pin->GetReal("problem", "T0");
  Tmin = pin->GetReal("problem", "Tmin");
  clat = pin->GetOrAddReal("problem", "clat", 0.);
  xH2S = pin->GetReal("problem", "xH2S");
  metallicity = pin->GetOrAddReal("problem", "metallicity", 0.);
  xNa = pin->GetReal("problem", "xNa");
  xNa *= pow(10., metallicity);
  xKCl = pin->GetReal("problem", "xKCl");
  xKCl *= pow(10., metallicity);
}

void RadiationBand::addAbsorber(std::string name, std::string file, ParameterInput *pin)
{
  std::stringstream msg;

  xHe = pin->GetReal("problem", "xHe");
  xCH4 = pin->GetReal("problem", "xCH4");

  if (name == "mw_CIA") {
    pabs->addAbsorber(MwrAbsorberCIA(this, xHe, xCH4));
  } else if (name == "mw_NH3") {
    pabs->addAbsorber(MwrAbsorberNH3(this, {iNH3, iH2O}, xHe).SetModelHanley());
  } else if (name == "mw_H2O") {
    pabs->addAbsorber(MwrAbsorberH2O(this, iH2O, xHe));
  } else if (name == "mw_electron") {
    pabs->addAbsorber(MwrAbsorberElectron(this, ion));
  } else {
    msg << "### FATAL ERROR in RadiationBand::addAbsorber"
        << std::endl << "unknow absorber: '" << name <<"' ";
    ATHENA_ERROR(msg);
  }
}

#ifdef REAL_GAS_CP
void update_gamma(Real& gamma, Real const q[]) {
  Real T = q[IDN], cp_h2, cp_he, cp_ch4;
  if (T < 300.)
    cp_h2 = Hydrogen::cp_norm(T);
  else
    cp_h2 = Hydrogen::cp_nist(T);
  cp_he = Helium::cp_nist(T);
  cp_ch4 = Methane::cp_nist(T);

  Real cp_real = (1. - xHe - xCH4)*cp_h2 + xHe*cp_he + xCH4*cp_ch4;
  gamma = cp_real/(cp_real - Thermodynamics::Rgas);
}
#endif

// @sect3{Initial condition}

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  static_assert(HYDROSTATIC, "This problem requires turning on hydrostatic option");
  pdebug->Enter("ProblemGenerator: juno_mwr");
  std::stringstream &msg = pdebug->msg;

  //ReadJunoMWRProfile("Juno_MWR_PJ1345689_m24-m16_avgcoeff.fits", coeff, cov);
  //ReadWriteGeminiTEXESMap("Gemini_TEXES_GRS_2017_product.dat", coeff, iNH3);

  // Construct a 1D pseudo-moist adiabat at (T0,P0)
  Real gamma = pin->GetReal("hydro", "gamma");
  Real Rd = pthermo->GetRd();
  Real cp = gamma/(gamma - 1.)*Rd;

  Real x1min = pmy_mesh->mesh_size.x1min;
  Real x1max = pmy_mesh->mesh_size.x1max;
  Real H0 = phydro->scale_height;

  // estimate surface temperature and pressure
  // coordinate range needs to be adjusted for hydrostatic model
  Real Ps = P0*exp(-x1min/H0);
  Real Ts = T0*pow(Ps/P0,Rd/cp);
  Real dlnp, **w1, *z1, *t1;
  dlnp = (x1max - x1min)/pmy_mesh->mesh_size.nx1/(2.*H0);
  int nx1 = (int)((x1max - x1min)/(H0*dlnp));
  NewCArray(w1, nx1, NHYDRO+2*NVAPOR);
  z1 = new Real [nx1];
  t1 = new Real [nx1];

  int max_iter = 200, iter = 0;

  for (int n = 1; n <= NVAPOR; ++n) {
    Real qv = pin->GetReal("problem", "qvapor" + std::to_string(n))/1.E3;
    w1[0][n] = qv;
  }
  z1[0] = x1min;
  for (int i = 1; i < nx1; ++i)
    z1[i] = z1[i-1] + H0*dlnp;

  Real t0;
  if (Globals::my_rank == 0)
    std::cout << "- request T = " << T0 << " at P = " << P0 << std::endl;
  while (iter++ < max_iter) {
    pthermo->ConstructAtmosphere(w1, Ts, Ps, grav, -dlnp, nx1, Adiabat::pseudo, 1.);

    // Replace adiabatic atmosphere with isothermal atmosphere if temperature is too low
    int ii = 0;
    for (; ii < nx1; ++ii)
      if (pthermo->GetTemp(w1[ii]) < Tmin) break;
    Real Tv = w1[ii][IPR]/(w1[ii][IDN]*Rd);
    for (int i = ii; i < nx1; ++i) {
      w1[i][IDN] = w1[i][IPR]/(Rd*Tv);
      for (int n = 1; n <= NVAPOR; ++n)
        w1[i][n] = w1[ii][n];
    }

    // Find T at p = p0
    for (int i = 0; i < nx1; ++i) {
      t1[i] = pthermo->GetTemp(w1[i]);
    }
    t0 = interp1(0, t1, z1, nx1);

    Ts += T0 - t0;
    if (Globals::my_rank == 0)
      std::cout << "- iteration #" << iter << ": " << "T = " << t0  << std::endl;
    if (fabs(T0 - t0) < 0.01) break;
  }

  if (iter > max_iter) {
    msg << "### FATAL ERROR in problem generator"
        << std::endl << "maximum iteration reached."
        << std::endl << "T0 = " << t0;
    ATHENA_ERROR(msg);
  }

  // Change to log quantity
  for (int i = 0; i < nx1; ++i)
    for (int n = 0; n < NHYDRO+2*NVAPOR; ++n) {
      if (n == IVX || n == IVY || n == IVZ)
        w1[i][n] = 0.;
      else
        w1[i][n] = log(w1[i][n]);
    }

  // setup initial condition
  // Z = -H0*log(P/P0)
  // log(P) = log(P0) - Z/H0
  for (int i = is; i <= ie; ++i) {
    Real buf[NHYDRO+2*NVAPOR];
    // set gas concentration and velocity
    interpn(buf, &pcoord->x1v(i), *w1, z1, &nx1, 1, NHYDRO+2*NVAPOR);
    for (int n = 0; n < NHYDRO; ++n)
      for (int k = ks; k <= ke; ++k)
        for (int j = js; j <= je; ++j) {
          if (n == IVX || n == IVY || n == IVZ)
            phydro->w(n,k,j,i) = 0.;
          else
            phydro->w(n,k,j,i) = exp(buf[n]);
        }

    // set chemical tracer, electron and Na
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j) {
        Real temp = pthermo->GetTemp(phydro->w.at(k,j,i));
        Real pH2S = xH2S*phydro->w(IPR,k,j,i);
        Real pNa = xNa*phydro->w(IPR,k,j,i);
        Real svp = sat_vapor_p_Na_H2S_Visscher(temp, pH2S);
        pNa = std::min(svp, pNa);
        pscalars->s(iNa,k,j,i) = pNa/(Thermodynamics::kBoltz*temp);
        pscalars->s(ion,k,j,i) = saha_ionization_electron_density(temp, pscalars->s(iNa,k,j,i), 5.14);
      }
  }

  // replace adiabatic temperature profile with measurement above 1 bar
  if (pin->GetOrAddBoolean("problem", "use_fletcher16_cirs", false)) {
    // clat->glat, pa -> bar
    Real glat = jup_centric2graphic(clat);
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i) {
          if (phydro->w(IPR,k,j,i) > 1.E5) continue;
          Real temp_ad = pthermo->GetTemp(phydro->w.at(k,j,i));
          Real temp_real = Jupiter::get_temp_fletcher16_cirs(glat, phydro->w(IPR,k,j,i)/1.E5);
          Real R = pthermo->RovRd(phydro->w.at(k,j,i))*Rd;
          phydro->w(IDN,k,j,i) = phydro->w(IPR,k,j,i)/(R*std::max(temp_real, temp_ad));
        }
  }

  // Microwave radiative transfer needs temperatures at cell interfaces, which are
  // interpolated from cell centered hydrodynamic variables.
  // Normally, the boundary conditions are taken care of internally.
  // But, since we call radiative tranfer directly in pgen, we would need to update the
  // boundary conditions manually. The following lines of code updates the boundary
  // conditions.
  phydro->hbvar.SwapHydroQuantity(phydro->w, HydroBoundaryQuantity::prim);
  pbval->ApplyPhysicalBoundaries(0., 0.);

  // read profile updates from input
  std::vector<Real> Tp, NH3p;

  Tp = Vectorize<Real>(pin->GetString("problem", "Tp").c_str());
  NH3p = Vectorize<Real>(pin->GetString("problem", "NH3p").c_str());
  int nsample = pfit->plevel.size() - 2;   // excluding boundary condition

  Real **XpSample;
  NewCArray(XpSample, 1+NVAPOR, nsample+2);
  std::fill(*XpSample, *XpSample + (1+NVAPOR)*(nsample+2), 0.);

  for (int i = 0; i < nsample; ++i) {
    XpSample[0][i+1] = Tp[i];
    XpSample[iNH3][i+1] = NH3p[i]*1.E-3;  // g/kg -> kg/kg
  }

  // Revise atmospheric profiles and store them in position (k,je).
  for (int k = ks; k <= ke; ++k)
    pfit->UpdateAtmosphere(&XpSample, k);

  // Swap the revised profile and baseline profile (k,je) <-> (k,js).
  for (int n = 0; n < NHYDRO; ++n)
    for (int k = ks; k <= ke; ++k)
      for (int i = is; i <= ie; ++i) {
        Real tmp = phydro->w(n,k,js,i);
        phydro->w(n,k,js,i) = phydro->w(n,k,je,i);
        phydro->w(n,k,je,i) = tmp;
      }

  // Calculate microwave radiation at time t = 0
  // These values maintain for all future time steps (t > 0,k,js)
  if (Globals::my_rank == 0)
    std::cout << "- running initial RT for all models" << std::endl;
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      prad->calculateRadiances(phydro->w, 0., k, j, is, ie+1);

  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);
  // Leave debug stack and print out all debug info.
  pdebug->Leave();

  // Clean up temporary memory
  FreeCArray(w1);
  FreeCArray(XpSample);
  delete[] z1;
  delete[] t1;
}
