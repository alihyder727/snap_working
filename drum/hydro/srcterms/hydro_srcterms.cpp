//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//  \brief Class to implement source terms in the hydro equations

// C headers

// C++ headers
#include <cstring>    // strcmp
#include <iostream>
#include <sstream>
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "../../parameter_input.hpp"
#include "../../debugger/debugger.hpp"
#include "../hydro.hpp"
#include "hydro_srcterms.hpp"

// HydroSourceTerms constructor

HydroSourceTerms::HydroSourceTerms(Hydro *phyd, ParameterInput *pin) {
  Debugger *pdbg =  phyd->pmy_block->pdebug;
  pdbg->Enter("HydroSouceTerms");
  //ATHENA_LOG("HydroSourceTerms");
  pmy_hydro_ = phyd;
  hydro_sourceterms_defined = false;
  std::stringstream &msg = pdbg->msg;

  // read point mass or constant acceleration parameters from input block

  // set the point source only when the coordinate is spherical or 2D
  // cylindrical.
  gm_ = pin->GetOrAddReal("problem","GM",0.0);
  if (gm_ != 0.0) {
    if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0
        || (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0
            && phyd->pmy_block->block_size.nx3==1)) {
      hydro_sourceterms_defined = true;
    } else {
      msg << "### FATAL ERROR in HydroSourceTerms constructor" << std::endl
          << "The point mass gravity works only in spherical polar coordinates"
          << "or in 2D cylindrical coordinates." << std::endl
          << "Check <problem> GM parameter in the input file." << std::endl;
      ATHENA_ERROR(msg);
    }
  }
  g1_ = pin->GetOrAddReal("hydro","grav_acc1",0.0);
  if (g1_ > 0.) ATHENA_WARN("grav_acc1 is positive.");
  if (g1_ != 0.0) {
    hydro_sourceterms_defined = true;
    msg << "- g1 = " << g1_ << std::endl;
  }

  g2_ = pin->GetOrAddReal("hydro","grav_acc2",0.0);
  if (g2_ != 0.0) {
    hydro_sourceterms_defined = true;
    msg << "- g2 = " << g1_ << std::endl;
  }

  g3_ = pin->GetOrAddReal("hydro","grav_acc3",0.0);
  if (g3_ != 0.0) {
    hydro_sourceterms_defined = true;
    msg << "- g3 = " << g1_ << std::endl;
  }

  // coriolis acceleration
  omega1_ = pin->GetOrAddReal("hydro","coriolis_acc1",0.0);
  if (omega1_ != 0.0) hydro_sourceterms_defined = true;

  omega2_ = pin->GetOrAddReal("hydro","coriolis_acc2",0.0);
  if (omega2_ != 0.0) hydro_sourceterms_defined = true;

  omega3_ = pin->GetOrAddReal("hydro","coriolis_acc3",0.0);
  if (omega3_ != 0.0) hydro_sourceterms_defined = true;

  omegax_ = pin->GetOrAddReal("hydro","OmegaX",0.0);
  if (omegax_ != 0.0) {
    hydro_sourceterms_defined = true;
    msg << "- Rotation in X-direction: " << omegax_ << std::endl;
  }

  omegay_ = pin->GetOrAddReal("hydro","OmegaY",0.0);
  if (omegay_ != 0.0) {
    hydro_sourceterms_defined = true;
    msg << "- Rotation in Y-direction: " << omegay_ << std::endl;
  }

  omegaz_ = pin->GetOrAddReal("hydro","OmegaZ",0.0);
  if (omegaz_ != 0.0) {
    hydro_sourceterms_defined = true;
    msg << "- Rotation in Z-direction: " << omegaz_ << std::endl;
  }

  // read shearing box parameters from input block
  Omega_0_ = pin->GetOrAddReal("problem","Omega0",0.0);
  qshear_  = pin->GetOrAddReal("problem","qshear",0.0);
  ShBoxCoord_ = pin->GetOrAddInteger("problem","shboxcoord",1);
  if ((Omega_0_ !=0.0) && (qshear_ != 0.0)) hydro_sourceterms_defined = true;

  if (SELF_GRAVITY_ENABLED) hydro_sourceterms_defined = true;

  UserSourceTerm = phyd->pmy_block->pmy_mesh->UserSourceTerm_;
  if (UserSourceTerm != nullptr) hydro_sourceterms_defined = true;
  pdbg->Leave();
}

//----------------------------------------------------------------------------------------
//! \fn void HydroSourceTerms::AddHydroSourceTerms
//  \brief Adds source terms to conserved variables

void HydroSourceTerms::AddHydroSourceTerms(const Real time, const Real dt,
                                           const AthenaArray<Real> *flux,
                                           const AthenaArray<Real> &prim,
                                           const AthenaArray<Real> &bcc,
                                           AthenaArray<Real> &cons) {
  MeshBlock *pmb = pmy_hydro_->pmy_block;

  // accleration due to point mass (MUST BE AT ORIGIN)
  if (gm_ != 0.0) PointMass(dt, flux, prim, cons);

  // constant acceleration (e.g. for RT instability)
  if (g1_ != 0.0 || g2_ != 0.0 || g3_ != 0.0) ConstantAcceleration(dt, flux,
                                                                   prim,cons);

  // coriolis acceleration in the axial direction
  if (omega1_ != 0.0 || omega2_ != 0.0 || omega3_ != 0.0)
    Coriolis123(dt, flux, prim, cons);

  // coriolis acceleration in the cartesian direction
  if (omegax_ != 0.0 || omegay_ != 0.0 || omegaz_ != 0.0)
    CoriolisXYZ(dt, flux, prim, cons);

  // Add new source terms here
  if (SELF_GRAVITY_ENABLED) SelfGravity(dt, flux, prim, cons);
  // shearing box source terms: tidal and Coriolis forces
  if ((Omega_0_ !=0.0) && (qshear_ != 0.0)) ShearingBoxSourceTerms(dt, flux,
                                                                   prim, cons);
  // MyNewSourceTerms()

  //  user-defined source terms
  if (UserSourceTerm != nullptr)
    UserSourceTerm(pmb, time,dt,prim,bcc,cons);

  return;
}
