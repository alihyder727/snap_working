//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mesh.cpp
//  \brief implementation of functions in MeshBlock class

// C headers

// C++ headers
#include <algorithm>  // sort()
#include <cstdlib>
#include <cstring>    // memcpy()
#include <ctime>      // clock(), CLOCKS_PER_SEC, clock_t
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../fft/athena_fft.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../gravity/gravity.hpp"
#include "../gravity/mg_gravity.hpp"
#include "../hydro/hydro.hpp"
#include "../parameter_input.hpp"
#include "../reconstruct/reconstruction.hpp"
#include "../scalars/scalars.hpp"
#include "../utils/buffer_utils.hpp"
#include "mesh.hpp"
#include "mesh_refinement.hpp"
#include "meshblock_tree.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "../chemistry/chemistry.hpp"
#include "../radiation/radiation.hpp"
#include "../physics/physics.hpp"
#include "../diagnostics/diagnostics.hpp"
#include "../debugger/debugger.hpp"
#include "../particles/particles.hpp"
#include "../inversion/inversion.hpp"
#include "../communicator/communicator.hpp"
#include "../turbulence/turbulence_model.hpp"

//----------------------------------------------------------------------------------------
// MeshBlock constructor: constructs coordinate, boundary condition, hydro, field
//                        and mesh refinement objects.

MeshBlock::MeshBlock(int igid, int ilid, LogicalLocation iloc, RegionSize input_block,
                     BoundaryFlag *input_bcs, Mesh *pm, ParameterInput *pin,
                     int igflag, bool ref_flag) :
    pmy_mesh(pm), loc(iloc), block_size(input_block),
    gid(igid), lid(ilid), gflag(igflag), nuser_out_var(), prev(nullptr), next(nullptr),
    new_block_dt_{}, new_block_dt_hyperbolic_{}, new_block_dt_parabolic_{},
    new_block_dt_user_{},
    nreal_user_meshblock_data_(), nint_user_meshblock_data_(), cost_(1.0) {
  // initialize grid indices
  is = NGHOST;
  ie = is + block_size.nx1 - 1;

  ncells1 = block_size.nx1 + 2*NGHOST;
  ncc1 = block_size.nx1/2 + 2*NGHOST;
  if (pmy_mesh->f2) {
    js = NGHOST;
    je = js + block_size.nx2 - 1;
    ncells2 = block_size.nx2 + 2*NGHOST;
    ncc2 = block_size.nx2/2 + 2*NGHOST;
  } else {
    js = je = 0;
    ncells2 = 1;
    ncc2 = 1;
  }

  if (pmy_mesh->f3) {
    ks = NGHOST;
    ke = ks + block_size.nx3 - 1;
    ncells3 = block_size.nx3 + 2*NGHOST;
    ncc3 = block_size.nx3/2 + 2*NGHOST;
  } else {
    ks = ke = 0;
    ncells3 = 1;
    ncc3 = 1;
  }

  if (pm->multilevel) {
    cnghost = (NGHOST + 1)/2 + 1;
    cis = NGHOST; cie = cis + block_size.nx1/2 - 1;
    cjs = cje = cks = cke = 0;
    if (pmy_mesh->f2) // 2D or 3D
      cjs = NGHOST, cje = cjs + block_size.nx2/2 - 1;
    if (pmy_mesh->f3) // 3D
      cks = NGHOST, cke = cks + block_size.nx3/2 - 1;
  }

  // (probably don't need to preallocate space for references in these vectors)
  vars_cc_.reserve(3);
  vars_fc_.reserve(3);

  // construct objects stored in MeshBlock class.  Note in particular that the initial
  // conditions for the simulation are set in problem generator called from main, not
  // in the Hydro constructor

  // mesh-related objects
  // debug
  pdebug = new Debugger(this);

  // Boundary
  pbval  = new BoundaryValues(this, input_bcs, pin);

  // Coordinates
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
    pcoord = new Cartesian(this, pin, false);
  } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    pcoord = new Cylindrical(this, pin, false);
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    pcoord = new SphericalPolar(this, pin, false);
  } else if (std::strcmp(COORDINATE_SYSTEM, "minkowski") == 0) {
    pcoord = new Minkowski(this, pin, false);
  } else if (std::strcmp(COORDINATE_SYSTEM, "schwarzschild") == 0) {
    pcoord = new Schwarzschild(this, pin, false);
  } else if (std::strcmp(COORDINATE_SYSTEM, "kerr-schild") == 0) {
    pcoord = new KerrSchild(this, pin, false);
  } else if (std::strcmp(COORDINATE_SYSTEM, "gr_user") == 0) {
    pcoord = new GRUser(this, pin, false);
  }

  // Reconstruction: constructor may implicitly depend on Coordinates, and PPM variable
  // floors depend on EOS, but EOS isn't needed in Reconstruction constructor-> this is ok
  precon = new Reconstruction(this, pin);

  if (pm->multilevel) pmr = new MeshRefinement(this, pin);

  // physics-related, per-MeshBlock objects: may depend on Coordinates for diffusion
  // terms, and may enroll quantities in AMR and BoundaryVariable objs. in BoundaryValues

  // TODO(felker): prepare this section of the MeshBlock ctor to become more complicated
  // for several extensions:
  // 1) allow solver to compile without a Hydro class (or with a Hydro class for the
  // background fluid that is not dynamically evolved)
  // 2) MPI ranks containing MeshBlocks that solve a subset of the physics, e.g. Gravity
  // but not Hydro.
  // 3) MAGNETIC_FIELDS_ENABLED, SELF_GRAVITY_ENABLED, NSCALARS, (future) FLUID_ENABLED,
  // etc. become runtime switches

  // if (FLUID_ENABLED) {
    // if (this->hydro_block)
    phydro = new Hydro(this, pin);
    // } else
    // }
    // Regardless, advance MeshBlock's local counter (initialized to bvars_next_phys_id=1)
    // Greedy reservation of phys IDs (only 1 of 2 needed for Hydro if multilevel==false)
    pbval->AdvanceCounterPhysID(HydroBoundaryVariable::max_phys_id);
    //  }
  if (MAGNETIC_FIELDS_ENABLED) {
    // if (this->field_block)
    pfield = new Field(this, pin);
    pbval->AdvanceCounterPhysID(FaceCenteredBoundaryVariable::max_phys_id);
  }
  if (SELF_GRAVITY_ENABLED) {
    // if (this->grav_block)
    pgrav = new Gravity(this, pin);
    pbval->AdvanceCounterPhysID(CellCenteredBoundaryVariable::max_phys_id);
    if (SELF_GRAVITY_ENABLED == 2)
      pmg = new MGGravity(pmy_mesh->pmgrd, this);
  }
  if (NSCALARS > 0) {
    // if (this->scalars_block)
    pscalars = new PassiveScalars(this, pin);
    pbval->AdvanceCounterPhysID(CellCenteredBoundaryVariable::max_phys_id);
  }
  // KGF: suboptimal solution, since developer must copy/paste BoundaryVariable derived
  // class type that is used in each PassiveScalars, Gravity, Field, Hydro, ... etc. class
  // in order to correctly advance the BoundaryValues::bvars_next_phys_id_ local counter.

  // TODO(felker): check that local counter pbval->bvars_next_phys_id_ agrees with shared
  // Mesh::next_phys_id_ counter (including non-BoundaryVariable / per-MeshBlock reserved
  // values). Compare both private member variables via BoundaryValues::CheckCounterPhysID

  peos = new EquationOfState(this, pin);
  pthermo = new Thermodynamics(this, pin);
  prad = new Radiation(this, pin);
  pphy = new Physics(this, pin);
  // initialize particles and remove the first one (header)
  ppart = new Particles(this, pin);
  if (ppart->next != nullptr) {
    ppart = ppart->next;
    delete ppart->prev;
  } else {
    delete ppart;
    ppart = nullptr;
  }
  pchem = new Chemistry(this, pin);
  pdiag = new Diagnostics(this, pin);
  pfit = nullptr;
  pcomm = new Communicator(this);
  if (std::strcmp(TURBULENCE_MODEL, "KEpsilon") == 0) {
    pturb = new KEpsilonTurbulence(this, pin);
    pbval->AdvanceCounterPhysID(CellCenteredBoundaryVariable::max_phys_id);
  } else {
    pturb = new TurbulenceModel(this, pin);
  }

  // Create user mesh data
  InitUserMeshBlockData(pin);

  return;
}

//----------------------------------------------------------------------------------------
// MeshBlock constructor for restarts

MeshBlock::MeshBlock(int igid, int ilid, Mesh *pm, ParameterInput *pin,
                     LogicalLocation iloc, RegionSize input_block,
                     BoundaryFlag *input_bcs,
                     double icost, char *mbdata, int igflag) :
    pmy_mesh(pm), loc(iloc), block_size(input_block),
    gid(igid), lid(ilid), gflag(igflag), nuser_out_var(), prev(nullptr), next(nullptr),
    new_block_dt_{}, new_block_dt_hyperbolic_{}, new_block_dt_parabolic_{},
    new_block_dt_user_{},
    nreal_user_meshblock_data_(), nint_user_meshblock_data_(), cost_(icost) {
  // initialize grid indices
  is = NGHOST;
  ie = is + block_size.nx1 - 1;

  ncells1 = block_size.nx1 + 2*NGHOST;
  ncc1 = block_size.nx1/2 + 2*NGHOST;
  if (pmy_mesh->f2) {
    js = NGHOST;
    je = js + block_size.nx2 - 1;
    ncells2 = block_size.nx2 + 2*NGHOST;
    ncc2 = block_size.nx2/2 + 2*NGHOST;
  } else {
    js = je = 0;
    ncells2 = 1;
    ncc2 = 1;
  }

  if (pmy_mesh->f3) {
    ks = NGHOST;
    ke = ks + block_size.nx3 - 1;
    ncells3 = block_size.nx3 + 2*NGHOST;
    ncc3 = block_size.nx3/2 + 2*NGHOST;
  } else {
    ks = ke = 0;
    ncells3 = 1;
    ncc3 = 1;
  }

  if (pm->multilevel) {
    cnghost = (NGHOST + 1)/2 + 1;
    cis = NGHOST; cie = cis + block_size.nx1/2 - 1;
    cjs = cje = cks = cke = 0;
    if (pmy_mesh->f2) // 2D or 3D
      cjs = NGHOST, cje = cjs + block_size.nx2/2 - 1;
    if (pmy_mesh->f3) // 3D
      cks = NGHOST, cke = cks + block_size.nx3/2 - 1;
  }

  // (re-)create mesh-related objects in MeshBlock
  // debug
  pdebug = new Debugger(this);

  // Boundary
  pbval = new BoundaryValues(this, input_bcs, pin);

  // Coordinates
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
    pcoord = new Cartesian(this, pin, false);
  } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    pcoord = new Cylindrical(this, pin, false);
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    pcoord = new SphericalPolar(this, pin, false);
  } else if (std::strcmp(COORDINATE_SYSTEM, "minkowski") == 0) {
    pcoord = new Minkowski(this, pin, false);
  } else if (std::strcmp(COORDINATE_SYSTEM, "schwarzschild") == 0) {
    pcoord = new Schwarzschild(this, pin, false);
  } else if (std::strcmp(COORDINATE_SYSTEM, "kerr-schild") == 0) {
    pcoord = new KerrSchild(this, pin, false);
  } else if (std::strcmp(COORDINATE_SYSTEM, "gr_user") == 0) {
    pcoord = new GRUser(this, pin, false);
  }

  // Reconstruction (constructor may implicitly depend on Coordinates)
  precon = new Reconstruction(this, pin);

  if (pm->multilevel) pmr = new MeshRefinement(this, pin);

  // (re-)create physics-related objects in MeshBlock

  // if (FLUID_ENABLED) {
  // if (this->hydro_block)
  phydro = new Hydro(this, pin);
  // } else
  // }
  // Regardless, advance MeshBlock's local counter (initialized to bvars_next_phys_id=1)
  // Greedy reservation of phys IDs (only 1 of 2 needed for Hydro if multilevel==false)
  pbval->AdvanceCounterPhysID(HydroBoundaryVariable::max_phys_id);
  //  }
  if (MAGNETIC_FIELDS_ENABLED) {
    // if (this->field_block)
    pfield = new Field(this, pin);
    pbval->AdvanceCounterPhysID(FaceCenteredBoundaryVariable::max_phys_id);
  }
  if (SELF_GRAVITY_ENABLED) {
    // if (this->grav_block)
    pgrav = new Gravity(this, pin);
    pbval->AdvanceCounterPhysID(CellCenteredBoundaryVariable::max_phys_id);
    if (SELF_GRAVITY_ENABLED == 2)
      pmg = new MGGravity(pmy_mesh->pmgrd, this);
  }

  if (NSCALARS > 0) {
    // if (this->scalars_block)
    pscalars = new PassiveScalars(this, pin);
    pbval->AdvanceCounterPhysID(CellCenteredBoundaryVariable::max_phys_id);
  }

  peos = new EquationOfState(this, pin);
  pthermo = new Thermodynamics(this, pin);
  prad = new Radiation(this, pin);
  pphy = new Physics(this, pin);
  // initialize particles and remove the first one (header)
  ppart = new Particles(this, pin);
  if (ppart->next != nullptr) {
    ppart = ppart->next;
    delete ppart->prev;
  } else {
    delete ppart;
    ppart = nullptr;
  }
  pchem = new Chemistry(this, pin);
  pdiag = new Diagnostics(this, pin);
  pfit = nullptr;
  pcomm = new Communicator(this);
  if (std::strcmp(TURBULENCE_MODEL, "KEpsilon") == 0) {
    pturb = new KEpsilonTurbulence(this, pin);
    pbval->AdvanceCounterPhysID(CellCenteredBoundaryVariable::max_phys_id);
  } else {
    pturb = new TurbulenceModel(this, pin);
  }

  InitUserMeshBlockData(pin);

  std::size_t os = 0;
  // NEW_OUTPUT_TYPES:

  // load hydro and field data
  std::memcpy(phydro->u.data(), &(mbdata[os]), phydro->u.GetSizeInBytes());
  // outflow boundary condition
  if (pbval->block_bcs[inner_x1] == BoundaryFlag::outflow) {
    int kl = block_size.nx3 == 1 ? ks : ks-NGHOST;
    int ku = block_size.nx3 == 1 ? ke : ke+NGHOST;
    int jl = block_size.nx2 == 1 ? js : js-NGHOST;
    int ju = block_size.nx2 == 1 ? je : je+NGHOST;
    peos->ConservedToPrimitive(phydro->u, phydro->w1, pfield->b,
                                          phydro->w, pfield->bcc, pcoord,
                                          is-NGHOST, is-1, jl, ju, kl, ku);
    phydro->w1 = phydro->w;
  }
  // load it into the other memory register(s) too
  std::memcpy(phydro->u1.data(), &(mbdata[os]), phydro->u1.GetSizeInBytes());
  os += phydro->u.GetSizeInBytes();
  if (GENERAL_RELATIVITY) {
    std::memcpy(phydro->w.data(), &(mbdata[os]), phydro->w.GetSizeInBytes());
    os += phydro->w.GetSizeInBytes();
    std::memcpy(phydro->w1.data(), &(mbdata[os]), phydro->w1.GetSizeInBytes());
    os += phydro->w1.GetSizeInBytes();
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    std::memcpy(pfield->b.x1f.data(), &(mbdata[os]), pfield->b.x1f.GetSizeInBytes());
    std::memcpy(pfield->b1.x1f.data(), &(mbdata[os]), pfield->b1.x1f.GetSizeInBytes());
    os += pfield->b.x1f.GetSizeInBytes();
    std::memcpy(pfield->b.x2f.data(), &(mbdata[os]), pfield->b.x2f.GetSizeInBytes());
    std::memcpy(pfield->b1.x2f.data(), &(mbdata[os]), pfield->b1.x2f.GetSizeInBytes());
    os += pfield->b.x2f.GetSizeInBytes();
    std::memcpy(pfield->b.x3f.data(), &(mbdata[os]), pfield->b.x3f.GetSizeInBytes());
    std::memcpy(pfield->b1.x3f.data(), &(mbdata[os]), pfield->b1.x3f.GetSizeInBytes());
    os += pfield->b.x3f.GetSizeInBytes();
  }

  // (conserved variable) Passive scalars:
  if (NSCALARS > 0) {
    std::memcpy(pscalars->s.data(), &(mbdata[os]), pscalars->s.GetSizeInBytes());
    // load it into the other memory register(s) too
    std::memcpy(pscalars->s1.data(), &(mbdata[os]), pscalars->s1.GetSizeInBytes());
    os += pscalars->s.GetSizeInBytes();
  }

  // load user MeshBlock data
  for (int n=0; n<nint_user_meshblock_data_; n++) {
    std::memcpy(iuser_meshblock_data[n].data(), &(mbdata[os]),
                iuser_meshblock_data[n].GetSizeInBytes());
    os += iuser_meshblock_data[n].GetSizeInBytes();
  }
  for (int n=0; n<nreal_user_meshblock_data_; n++) {
    std::memcpy(ruser_meshblock_data[n].data(), &(mbdata[os]),
                ruser_meshblock_data[n].GetSizeInBytes());
    os += ruser_meshblock_data[n].GetSizeInBytes();
  }

  // load radiation data
  os += prad->loadRestartData(&(mbdata[os]));

  // load physics data
  os += pphy->LoadRestartData(&(mbdata[os]));

  // load particle data
  if (ppart != nullptr)
    os += ppart->LoadRestartData(&(mbdata[os]));

  // load turbulence data
  os += pturb->loadRestartData(&(mbdata[os]));

  return;
}

//----------------------------------------------------------------------------------------
// MeshBlock destructor

MeshBlock::~MeshBlock() {
  if (prev != nullptr) prev->next = next;
  if (next != nullptr) next->prev = prev;

  delete pcoord;
  delete precon;
  if (pmy_mesh->multilevel) delete pmr;

  delete phydro;
  if (MAGNETIC_FIELDS_ENABLED) delete pfield;
  delete peos;
  if (SELF_GRAVITY_ENABLED) delete pgrav;
  if (NSCALARS > 0) delete pscalars;

  // BoundaryValues should be destructed AFTER all BoundaryVariable objects are destroyed
  delete pbval;
  // delete user output variables array
  if (nuser_out_var > 0) {
    delete [] user_out_var_names_;
    delete [] user_out_var_longnames_;
    delete [] user_out_var_units_;
  }
  // delete user MeshBlock data
  if (nreal_user_meshblock_data_ > 0) delete [] ruser_meshblock_data;
  if (nint_user_meshblock_data_ > 0) delete [] iuser_meshblock_data;

  delete pthermo;
  delete prad;
  delete pphy;

  while (pdiag->prev != nullptr)
    delete pdiag->prev;
  while (pdiag->next != nullptr)
    delete pdiag->next;
  delete pdiag;

  while (pdebug->prev != nullptr)
    delete pdebug->prev;
  while (pdebug->next != nullptr)
    delete pdebug->next;
  delete pdebug;

  if (ppart != nullptr) {
    while (ppart->prev != nullptr)
      delete ppart->prev;
    while (ppart->next != nullptr)
      delete ppart->next;
    delete ppart;
  }
  delete pchem;

  if (pfit != nullptr)
    delete pfit;
  delete pcomm;
  delete pturb;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBlock::AllocateRealUserMeshBlockDataField(int n)
//  \brief Allocate Real AthenaArrays for user-defned data in MeshBlock

void MeshBlock::AllocateRealUserMeshBlockDataField(int n) {
  if (nreal_user_meshblock_data_ != 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in MeshBlock::AllocateRealUserMeshBlockDataField"
        << std::endl << "User MeshBlock data arrays are already allocated" << std::endl;
    ATHENA_ERROR(msg);
  }
  nreal_user_meshblock_data_ = n;
  ruser_meshblock_data = new AthenaArray<Real>[n];
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBlock::AllocateIntUserMeshBlockDataField(int n)
//  \brief Allocate integer AthenaArrays for user-defned data in MeshBlock

void MeshBlock::AllocateIntUserMeshBlockDataField(int n) {
  if (nint_user_meshblock_data_ != 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in MeshBlock::AllocateIntusermeshblockDataField"
        << std::endl << "User MeshBlock data arrays are already allocated" << std::endl;
    ATHENA_ERROR(msg);
    return;
  }
  nint_user_meshblock_data_=n;
  iuser_meshblock_data = new AthenaArray<int>[n];
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBlock::AllocateUserOutputVariables(int n)
//  \brief Allocate user-defined output variables

void MeshBlock::AllocateUserOutputVariables(int n) {
  if (n <= 0) return;
  if (nuser_out_var != 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in MeshBlock::AllocateUserOutputVariables"
        << std::endl << "User output variables are already allocated." << std::endl;
    ATHENA_ERROR(msg);
    return;
  }
  nuser_out_var = n;
  user_out_var.NewAthenaArray(nuser_out_var, ncells3, ncells2, ncells1);
  user_out_var_names_ = new std::string[n];
  user_out_var_longnames_ = new std::string[n];
  user_out_var_units_ = new std::string[n];
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MeshBlock::SetUserOutputVariableName(int n, const char *name)
//  \brief set the user-defined output variable name

void MeshBlock::SetUserOutputVariableName(int n, const char *name,
  const char *long_name, const char *units) {
  if (n >= nuser_out_var) {
    std::stringstream msg;
    msg << "### FATAL ERROR in MeshBlock::SetUserOutputVariableName"
        << std::endl << "User output variable is not allocated." << std::endl;
    ATHENA_ERROR(msg);
    return;
  }
  user_out_var_names_[n] = name;
  user_out_var_longnames_[n] = long_name;
  user_out_var_units_[n] = units;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn std::size_t MeshBlock::GetBlockSizeInBytes()
//  \brief Calculate the block data size required for restart.

std::size_t MeshBlock::GetBlockSizeInBytes() {
  std::size_t size;
  // NEW_OUTPUT_TYPES:
  size = phydro->u.GetSizeInBytes();
  if (GENERAL_RELATIVITY) {
    size += phydro->w.GetSizeInBytes();
    size += phydro->w1.GetSizeInBytes();
  }
  if (MAGNETIC_FIELDS_ENABLED)
    size += (pfield->b.x1f.GetSizeInBytes() + pfield->b.x2f.GetSizeInBytes()
             + pfield->b.x3f.GetSizeInBytes());
  if (SELF_GRAVITY_ENABLED)
    size += pgrav->phi.GetSizeInBytes();
  if (NSCALARS > 0)
    size += pscalars->s.GetSizeInBytes();

  // calculate user MeshBlock data size
  for (int n=0; n<nint_user_meshblock_data_; n++)
    size += iuser_meshblock_data[n].GetSizeInBytes();
  for (int n=0; n<nreal_user_meshblock_data_; n++)
    size += ruser_meshblock_data[n].GetSizeInBytes();

  // radiation data
  size += prad->getRestartDataSizeInBytes();

  // physics data
  size += pphy->RestartDataSizeInBytes();

  // particle data
  if (ppart != nullptr)
    size += ppart->RestartDataSizeInBytes();

  // turbulence data
  size += pturb->getRestartDataSizeInBytes();

  return size;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBlock::SetCostForLoadBalancing(double cost)
//  \brief stop time measurement and accumulate it in the MeshBlock cost

void MeshBlock::SetCostForLoadBalancing(double cost) {
  if (pmy_mesh->lb_manual_) {
    cost_ = std::min(cost, TINY_NUMBER);
    pmy_mesh->lb_flag_ = true;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBlock::ResetTimeMeasurement()
//  \brief reset the MeshBlock cost for automatic load balancing

void MeshBlock::ResetTimeMeasurement() {
  if (pmy_mesh->lb_automatic_) cost_ = TINY_NUMBER;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBlock::StartTimeMeasurement()
//  \brief start time measurement for automatic load balancing

void MeshBlock::StartTimeMeasurement() {
  if (pmy_mesh->lb_automatic_) {
#ifdef OPENMP_PARALLEL
    lb_time_ = omp_get_wtime();
#else
    lb_time_ = static_cast<double>(clock());
#endif
  }
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBlock::StartTimeMeasurement()
//  \brief stop time measurement and accumulate it in the MeshBlock cost

void MeshBlock::StopTimeMeasurement() {
  if (pmy_mesh->lb_automatic_) {
#ifdef OPENMP_PARALLEL
    lb_time_ = omp_get_wtime() - lb_time_;
#else
    lb_time_ = static_cast<double>(clock()) - lb_time_;
#endif
    cost_ += lb_time_;
  }
}


void MeshBlock::RegisterMeshBlockData(AthenaArray<Real> &pvar_cc) {
  vars_cc_.push_back(pvar_cc);
  return;
}


void MeshBlock::RegisterMeshBlockData(FaceField &pvar_fc) {
  vars_fc_.push_back(pvar_fc);
  return;
}

int MeshBlock::GetNumVariablesInCell() {
  int nvar = NHYDRO + NSCALARS;
  Particles *p = ppart;
  while (p != nullptr) {
    nvar += p->u.GetDim4();
    p = p->next;
  }
  return nvar;
}


// TODO(felker): consider merging the MeshRefinement::pvars_cc/fc_ into the
// MeshBlock::pvars_cc/fc_. Would need to weaken the MeshBlock std::vector to use tuples
// of pointers instead of a std::vector of references, so that:
// - nullptr can be passed for the second entry if multilevel==false
// - we can rebind the pointers to Hydro for GR purposes in bvals_refine.cpp
// If GR, etc. in the future requires additional flexiblity from non-refinement load
// balancing, we will need to use ptrs instead of references anyways, and add:

// void MeshBlock::SetHydroData(HydroBoundaryQuantity hydro_type)
//   Hydro *ph = pmy_block_->phydro;
//   // hard-coded assumption that, if multilevel, then Hydro is always present
//   // and enrolled in mesh refinement in the first pvars_cc_ vector entry
//   switch (hydro_type) {
//     case (HydroBoundaryQuantity::cons): {
//       pvars_cc_.front() = &ph->u;
//       break;
//     }
//     case (HydroBoundaryQuantity::prim): {
//       pvars_cc_.front() = &ph->w;
//       break;
//     }
//   }
//   return;
// }
