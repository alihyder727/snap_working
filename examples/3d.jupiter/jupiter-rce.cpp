// C/C++ header
#include <ctime>
#include <sstream>
#include <math.h>
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
#include "../reconstruct/interpolation.hpp"
// include variables: Real heating & Real heat_pres [ref: Guerlet et al. 2020]
#include "./heating_profile.txt"
#include "./primitive_1d.txt"
// MPI headers
#ifdef MPI_PARALLEL
  #include <mpi.h>
#endif

// molecules
enum {iH2O = 1, iNH3 = 2, iH2Oc = 3, iNH3c = 4, iH2Op = 5, iNH3p = 6};
Real grav, P0, T0, Tmin, prad, sponge_tau, sponge_width;
// internal heating and radiative heating parameters
// hflux: internal heat flux
Real hflux;
// forcing parameters
bool interp_pres;
// flags for horozontal differential cooling and vertical differential heating
bool VertDiffHeat, HoriDiffHeat, NewtonCoolFlag;
Real qH2O, qNH3;
Real Omega, center_deg, sponge_lat, NewtonCoolRate, 
     pdrag, Kv, beta_f, p_cool_bot, p_cool_top, 
     EqNoDragLat,qRelaxT;
int  Ndata;
void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{

  if (interp_pres) {
    Ndata = 7;
    AllocateUserOutputVariables(Ndata+NHYDRO);
    SetUserOutputVariableName(0,  "temp");
    SetUserOutputVariableName(1,  "theta");
    SetUserOutputVariableName(2,  "thetav");
    SetUserOutputVariableName(3,  "mse");
    SetUserOutputVariableName(4,  "rh_h2o");
    SetUserOutputVariableName(5,  "rh_nh3");
    SetUserOutputVariableName(6,  "phi");
    // NHYDRO-1+NVAPOR*NPHASE (pressure excluded)
    SetUserOutputVariableName(7,  "rho");
    SetUserOutputVariableName(8,  "vapor1");
    SetUserOutputVariableName(9,  "vapor2");
    SetUserOutputVariableName(10, "cloud11");
    SetUserOutputVariableName(11, "cloud21");
    SetUserOutputVariableName(12, "cloud12");
    SetUserOutputVariableName(13, "cloud22");
    SetUserOutputVariableName(14, "vel1");
    SetUserOutputVariableName(15, "vel2");
    SetUserOutputVariableName(16, "vel3");
    // pressure coordinate
    SetUserOutputVariableName(17, "pres_coord");
  } else {
    // height coordinate output
    AllocateUserOutputVariables(6);
    SetUserOutputVariableName(0, "temp");
    SetUserOutputVariableName(1, "theta");
    SetUserOutputVariableName(2, "thetav");
    SetUserOutputVariableName(3, "mse");
    SetUserOutputVariableName(4, "rh_h2o");
    SetUserOutputVariableName(5, "rh_nh3");
  }

  // save the 1D-pressure profile for interpolation and pressure coordinate output
  // constants
  int nx1 = pmy_mesh->mesh_size.nx1;
  int nx2 = pmy_mesh->mesh_size.nx2;
  int nx3 = pmy_mesh->mesh_size.nx3;
  // allocate array
  AllocateRealUserMeshBlockDataField(2);
  // initial pressure field
  ruser_meshblock_data[0].NewAthenaArray(ncells3,ncells2,ncells1);
  // initial temperature field
  ruser_meshblock_data[1].NewAthenaArray(ncells3,ncells2,ncells1);
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
{
  Real dq[1+NVAPOR], rh;

  // constants
  int nx1 = ie-is+1;  
  AthenaArray<Real> &x1flux = phydro->flux[X1DIR];
  AthenaArray<Real> &x2flux = phydro->flux[X2DIR];
  AthenaArray<Real> &x3flux = phydro->flux[X3DIR];

  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        user_out_var(0,k,j,i) = pthermo->Temp(phydro->w.at(k,j,i));
        user_out_var(1,k,j,i) = pthermo->Theta(phydro->w.at(k,j,i), P0);
        // theta_v
        user_out_var(2,k,j,i) = user_out_var(1,k,j,i)*pthermo->Qeps(phydro->w.at(k,j,i));
        user_out_var(3,k,j,i) = pthermo->MSE(phydro->w.at(k,j,i), grav*pcoord->x1v(i));
        pthermo->SaturationSurplus(dq, phydro->w.at(k,j,i));
        rh = phydro->w(iH2O,k,j,i)/(phydro->w(iH2O,k,j,i) - dq[iH2O]);
        user_out_var(4,k,j,i) = rh;
        rh = phydro->w(iNH3,k,j,i)/(phydro->w(iNH3,k,j,i) - dq[iNH3]);
        user_out_var(5,k,j,i) = rh;
      }
    }
  }
  // interpolate everything to the isobaric coordinate

  // alloacte 1d columes.
  Real *pre_1d, *data_1d;
  // center,left, and right values for Geopotential height
  AthenaArray<Real> Z,Zl,Zr;
  Z.NewAthenaArray(ncells3,ncells2,ncells1);
  Zl.NewAthenaArray(ncells3+1,ncells2+1,ncells1);
  Zr.NewAthenaArray(ncells3+1,ncells2+1,ncells1);
  data_1d   = new Real [nx1]; // 1d data array (under the height coordinate)
  pre_1d    = new Real [nx1]; // pressure colume (under the height coordinate)

  if (interp_pres) {
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = 0; i < nx1; ++i) {
        Real x1 = pcoord->x1v(i+is);
        user_out_var(6,k,j,i+is) = x1;
        // pressure colume (under the height coordinate)
        pre_1d[i] = phydro->w(IPR,k,j,i+is);
        // isobaric coordinate
        user_out_var(Ndata+NHYDRO-1,k,j,i+is) = ruser_meshblock_data[0](k,j,i+is);
      }
      // [1] thermodynamic quantaties: T,theta,thetav,MSE,rh1,rh2,phi
      for (int idata = 0; idata < Ndata; ++idata) {
        for (int i = 0; i < nx1; ++i) {
          // 1d data array (under the height coordinate)
          data_1d[i] = user_out_var(idata,k,j,i+is);
        } 
        for (int i = 0; i < nx1; ++i) {
          // data interpolation (isobaric coordinate)
          user_out_var(idata,k,j,i+is) = interp1(user_out_var(Ndata+NHYDRO-1,k,j,i+is),data_1d,pre_1d,nx1);
        }
      }
      // copy geopotential height from user_out_var to Z
      for (int i = is; i <= ie; ++i) {
        Z(k,j,i) = user_out_var(Ndata-1,k,j,i);
      }
      // [2] hydros: density, velocity, and mixing ratio of tracers
      for (int idata = Ndata; idata < Ndata+NHYDRO-1; ++idata) {
        for (int i = 0; i < nx1; ++i) {
          // 1d data array (under the height coordinate)
          data_1d[i] = phydro->w(idata-Ndata,k,j,i+is);
        }
        for (int i = 0; i < nx1; ++i) {
          // data interpolation (isobaric coordinate)
          user_out_var(idata,k,j,i+is) = interp1(user_out_var(Ndata+NHYDRO-1,k,j,i+is),data_1d,pre_1d,nx1);
        }
      }
    }
  }
  // [3] adjust u and v with geopotential height change under pressure coordinate.
  //     u_p = dx'/dt = (dx'/dx)*(dx/dt) = (dx'/dx)*u_z = (\sqrt(dZ^2+dx^2)/dx)*u_z
  int u_idx = Ndata + IVX;
  int v_idx = Ndata + IVY;
  // 1. zonal velocity (x3 direction)
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        // interpolate geopotential height
        Zl(k+1,j,i) = interp_weno5(Z(k+2,j,i),Z(k+1,j,i),Z(k,j,i),
                                   Z(k-1,j,i),Z(k-2,j,i)); // left
        Zr(k,j,i)   = interp_weno5(Z(k-2,j,i),Z(k-1,j,i),Z(k,j,i),
                                   Z(k+1,j,i),Z(k+2,j,i)); // right
        Real dx3 = pcoord->dx3v(k);
        Real Z_grad = Zr(k,j,i) - Zl(k+1,j,i);
        user_out_var(u_idx,k,j,i) *= sqrt(Z_grad*Z_grad+dx3*dx3)/dx3;
      }
    }
  }
  // 2. meridional velocity (x2 direction)
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        // interpolate geopotential height
        Zl(k,j+1,i) = interp_weno5(Z(k,j+2,i),Z(k,j+1,i),Z(k,j,i),
                                   Z(k,j-1,i),Z(k,j-2,i)); // left
        Zr(k,j,i)   = interp_weno5(Z(k,j-2,i),Z(k,j-1,i),Z(k,j,i),
                                   Z(k,j+1,i),Z(k,j+2,i)); // right
        Real dx2 = pcoord->dx2v(j);
        Real Z_grad = Zr(k,j,i) - Zl(k,j+1,i);
        user_out_var(v_idx,k,j,i) *= sqrt(Z_grad*Z_grad+dx2*dx2)/dx2;
      }
    }
  }}
  Z.DeleteAthenaArray();
  Zl.DeleteAthenaArray();
  Zr.DeleteAthenaArray();
  delete[] data_1d;
  delete[] pre_1d;
}

void Forcing(MeshBlock *pmb, Real const time, Real const dt,
    AthenaArray<Real> const &w, AthenaArray<Real> const &bcc, AthenaArray<Real> &u)
{
  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;
  Thermodynamics *pthermo = pmb->pthermo;

  // estimate the size of the domain
  Real dz = pmb->pcoord->x1v(is+1)-pmb->pcoord->x1v(is);

  Real x1max = pmb->pmy_mesh->mesh_size.x1max;
  Real x3min = pmb->pmy_mesh->mesh_size.x3min;
  Real x3max = pmb->pmy_mesh->mesh_size.x3max;
  Real Lx = x3max - x3min;

  Real x2min = pmb->pmy_mesh->mesh_size.x2min;
  Real x2max = pmb->pmy_mesh->mesh_size.x2max;
  Real Ly = x2max - x2min;

  Real MeshHeatWork;
  Real MeshHeatCapa;
  Real DomainHeatWork;
  Real DomainHeatCapa;
  Real hrate, LatHeatCoeff;

  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      if (pmb->pbval->apply_bndry_fn_[BoundaryFace::inner_x1]){
        // internal heat flux 
        Real HeatRate = hflux/dz; // [W/m^3]
        // only for the bottom cell
        u(IEN,k,j,is) += dt*HeatRate*(1.+
                         1.E-4*sin(2.*M_PI*rand()/RAND_MAX));
        // relax the vapor mixing ratio of the bottom cell to the initial value
        u(iH2O,k,j,is) += w(IDN,k,j,is)*(qH2O-w(iH2O,k,j,is))*dt*qRelaxT;
        u(iNH3,k,j,is) += w(IDN,k,j,is)*(qNH3-w(iNH3,k,j,is))*dt*qRelaxT;
      }
      // loop over vertical
      for (int i = is; i <= ie; ++i) {
        Real dx2 = pmb->pcoord->dx2v(j);
        Real dx3 = pmb->pcoord->dx3v(k);
        // Coriolis force
        Real x2 = pmb->pcoord->x2v(j);

        // artificial heating in the stratosphere
        Real temp = pthermo->Temp(w.at(k,j,i));
        if ( x1max-pmb->pcoord->x1v(i) < sponge_width ) {
          Real strato_heating = w(IDN,k,j,i)*cv*(Tmin - temp)*sponge_tau*dt;
          u(IEN,k,j,i) += strato_heating;
          MeshHeatWork += strato_heating*dz*dx2*dx3;
        }

        // top sponge layer
        if ( x1max-pmb->pcoord->x1v(i) < sponge_width ) {
          Real dw = w(IVX,k,j,i)*sponge_tau*dt;
          Real dv = w(IVY,k,j,i)*sponge_tau*dt;
          Real du = w(IVZ,k,j,i)*sponge_tau*dt;
          u(IM1,k,j,i) -= w(IDN,k,j,i)*dw;
        }

        // * insolation near the tropopause *
        // distance to the equator where f = 0
        Real dy2Eq = local_f/beta_f;
        // equivalent planetary radius
        Real eq_Rp = 2.*Omega/beta_f;
        // angle of insolation
        if (HoriDiffHeat) { 
          LatHeatCoeff = 1.-fabs(dy2Eq/eq_Rp);        
        } else {
          LatHeatCoeff = 1.;
        }
        // interpolate the heating profile
        if (VertDiffHeat) {
          hrate = interp1(w(IPR,k,j,i)/1.E5,heating,heat_pres,13); // [K/Jupiter day]
          hrate /= 35733.;
        } else {
          hrate = 0.;
        }
        if ((w(IPR,k,j,i) < p_cool_bot) && (x1max-pmb->pcoord->x1v(i) > sponge_width)){
          // total heat capacity in this meshblock
          MeshHeatCapa += w(IDN,k,j,i)*cv*dx2*dx3; // [J/(K*m)]
        }
        if (w(IPR,k,j,i) < 3.E5) {
          // total radiative heating in this meshblock
          MeshHeatWork += hrate*w(IDN,k,j,i)*cv*LatHeatCoeff*dz*dx2*dx3; // [W]
          // radiative heating
          u(IEN,k,j,i) += dt*hrate*w(IDN,k,j,i)*cv*LatHeatCoeff; // [J/m^3]
        }
      }
    }
  }

  // sum all body cooling across meshblocks: DomainCoolWork [W]
  #ifdef MPI_PARALLEL
    MPI_Allreduce(&MeshHeatWork,&DomainHeatWork,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&MeshHeatCapa,&DomainHeatCapa,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  #else
    DomainHeatWork = MeshHeatWork;
    DomainHeatCapa = MeshHeatCapa;
  #endif
  // estimate spatially constant radiative cooling flux 
  Real CoolFlux = -DomainHeatWork/(Lx*Ly*dt)-hflux; // [W/m^2]
  // debug output
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //if (rank == 0) { 
  //  std::cout << "Cooling Flux: " << CoolFlux 
  //            << ", Heating Flux:" << DomainHeatWork/(Lx*Ly*dt) << std::endl;
  //}

  Real AvgHeatCapa = DomainHeatCapa/(Lx*Ly);     // [J/(K*m^3)]
  Real crate = CoolFlux/(dz*AvgHeatCapa);        // [K/s]
  for (int k = ks; k <= ke; ++k) 
    for (int j = js; j <= je; ++j) 
      for (int i = is; i <= ie; ++i) {
        if ((w(IPR,k,j,i) < p_cool_bot) && (x1max-pmb->pcoord->x1v(i) > sponge_width)){
          Real cv = pthermo->Cv(w.at(k,j,i));
          if (NewtonCoolFlag) {
            // newtonian cooling
            u(IEN,k,j,i) -= dt*w(IDN,k,j,i)*cv*NewtonCoolRate*
                       (pthermo->Temp(w.at(k,j,i))-pmb->ruser_meshblock_data[1](k,j,i))*
                       (1.+1.E-4*sin(2.*M_PI*rand()/RAND_MAX));
          } else { 
            // constant radiative body cooling everywhere
            u(IEN,k,j,i) += dt*w(IDN,k,j,i)*cv*crate;
          }
        }
      }
 
}

Real TotalH2O(MeshBlock *pmb, int iout)
{
  Real tot_h2o = 0.;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;

  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
      for(int i=is; i<=ie; i++) {
        // total water mixing ratio * density * cell volume
        tot_h2o += (pmb->phydro->w(iH2O,k,j,i)+pmb->phydro->w(iH2Oc,k,j,i)+pmb->phydro->w(iH2Op,k,j,i))*
                    pmb->phydro->w(IDN,k,j,i)*pmb->pcoord->GetCellVolume(k,j,i);
      }
    }
  }

  return tot_h2o;
}

Real TotalNH3(MeshBlock *pmb, int iout)
{
  Real tot_nh3 = 0.;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;

  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
      for(int i=is; i<=ie; i++) {
        // total nh3 mixing ratio * density * cell volume
        tot_nh3 += (pmb->phydro->w(iNH3,k,j,i)+pmb->phydro->w(iNH3c,k,j,i)+pmb->phydro->w(iNH3p,k,j,i))*
                    pmb->phydro->w(IDN,k,j,i)*pmb->pcoord->GetCellVolume(k,j,i);
      }
    }
  }

  return tot_nh3;
}

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  grav = - pin->GetReal("hydro", "grav_acc1");
  P0 = pin->GetReal("problem", "P0");
  T0 = pin->GetReal("problem", "T0");
  Tmin = pin->GetReal("problem", "Tmin");
  prad = pin->GetReal("problem", "prad");
  hflux = pin->GetReal("problem", "hflux");  // internal heat flux [W/m^2]
  p_cool_bot = pin->GetReal("problem", "p_cool_bot");
  p_cool_top = pin->GetReal("problem", "p_cool_top");
  sponge_tau = pin->GetReal("problem", "sponge_tau");
  sponge_width = pin->GetReal("problem", "sponge_width");
  Omega = pin->GetReal("problem", "Omega");
  pdrag = pin->GetReal("problem", "pdrag");
  Kv = pin->GetReal("problem", "Kv")/86400;
  qRelaxT = pin->GetReal("problem", "qRelaxT");
  center_deg = pin->GetReal("problem", "center_deg");
  beta_f = pin->GetReal("problem", "beta_f");
  interp_pres = pin->GetOrAddBoolean("problem", "interp_pres", false);
  VertDiffHeat = pin->GetOrAddBoolean("problem", "VertDiffHeat", false);
  HoriDiffHeat = pin->GetOrAddBoolean("problem", "HoriDiffHeat", false);
  NewtonCoolFlag = pin->GetOrAddBoolean("problem", "NewtonCoolFlag", false);
  NewtonCoolRate = pin->GetReal("problem", "NewtonCoolRate");
  EqNoDragLat = pin->GetReal("problem", "EqNoDragLat");
  qH2O = pin->GetReal("problem", "qH2O");
  qNH3 = pin->GetReal("problem", "qNH3");

  EnrollUserExplicitSourceFunction(Forcing);
  AllocateUserHistoryOutput(2);
  EnrollUserHistoryOutput(0, TotalH2O, "tot_H2O");
  EnrollUserHistoryOutput(1, TotalNH3, "tot_NH3");
}

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  Real qH2O = pin->GetReal("problem", "qH2O");
  Real qNH3 = pin->GetReal("problem", "qNH3");
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
  NewCArray(w1, nx1, NHYDRO);
  z1 = new Real [nx1];
  p1 = new Real [nx1];
  t1 = new Real [nx1];

  // 1.1 estimate surface temperature and pressure
  Real Rd = pthermo->GetRd();
  Real cp = gamma/(gamma - 1.)*Rd;
  Real Ts = T0 - grav/cp*x1min;
  Real Ps = P0*pow(Ts/T0, cp/Rd);
  int max_iter = 200, iter = 0;

  w1[0][iH2O] = qH2O;
  w1[0][iNH3] = qNH3;
  z1[0] = x1min;
  for (int i = 1; i < nx1; ++i)
    z1[i] = z1[i-1] + dz;

  Real t0, p0;
  while (iter++ < max_iter) {
    pthermo->ConstructAdiabat(w1, Ts, Ps, grav, dz, nx1, Adiabat::pseudo);

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
    std::stringstream msg;
    msg << "### FATAL ERROR in problem generator"
        << std::endl << "maximum iteration reached."
        << std::endl << "T0 = " << t0
        << std::endl << "P0 = " << p0;
    ATHENA_ERROR(msg);
  }

  // setup initial condition
  int kl = block_size.nx3 == 1 ? ks : ks-NGHOST;
  int ku = block_size.nx3 == 1 ? ke : ke+NGHOST;
  int jl = block_size.nx2 == 1 ? js : js-NGHOST;
  int ju = block_size.nx2 == 1 ? je : je+NGHOST;
  for (int i = is; i <= ie; ++i) {
    Real buf[NHYDRO];
    interpn(buf, &pcoord->x1v(i), *w1, z1, &nx1, 1, NHYDRO);
    buf[IVX] = buf[IVY] = buf[IVZ] = 0.;
    // set cloud to zero
    for (int n = 1+NVAPOR; n < NMASS; ++n)
      buf[n] = 0.;
    for (int n = 0; n < NHYDRO; ++n)
      for (int k = kl; k <= ku; ++k)
        for (int j = jl; j <= ju; ++j)
          phydro->w(n,k,j,i) = buf[n];
  }

  // setup open lower boundary
  if (pbval->block_bcs[inner_x1] == BoundaryFlag::outflow) {
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j) {
        for (int n = 0; n < NHYDRO; ++n)
          w1[0][n] = phydro->w(n,k,j,is);
        Real P1 = w1[0][IPR];
        Real T1 = pthermo->Temp(w1[0]);
        Real dz = pcoord->dx1f(is);
        // adiabatic extrapolate half a grid down
        pthermo->ConstructAdiabat(w1, T1, P1, grav, -dz/2., 2, Adiabat::reversible);
        for (int n = 0; n < NHYDRO; ++n)
          for (int i = 1; i <= NGHOST; ++i)
            phydro->w(n,k,j,is-i) = w1[1][n];
      }
  }

  for(int k = ks; k <= ke; k++) {
    for(int j = js; j <= je; j++) {
      for (int i = is; i <= ie; i++) {
        // initialization
        // pressure
        ruser_meshblock_data[0](k,j,i) = phydro->w(IPR,k,j,i);
        // temperature
        ruser_meshblock_data[1](k,j,i) = pthermo->Temp(phydro->w.at(k,j,i));
      }
    }
  }

  srand(Globals::my_rank + time(0));
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);
  pphy->Initialize(phydro->w);

  FreeCArray(w1);
  delete[] z1;
  delete[] p1;
  delete[] t1;
}


