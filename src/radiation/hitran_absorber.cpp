// C/C++ headers
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <cstring>

// Athena++ headers
#include "../math/interpolation.h"
#include "hitran_absorber.hpp"

// External library headers
#ifdef NETCDFOUTPUT
  #include <netcdf.h>
#endif
std::ostream& operator<<(std::ostream &os, HitranAbsorber const& ab)
{
  for (int i = 0; i < 3; ++i) {
    os << "Axis " << i << " = " << ab.len_[i] << std::endl;
    if (i == 0) {
      os << "Minimum value = " << *std::min_element(ab.axis_.begin(), ab.axis_.begin() + ab.len_[i]) << std::endl;
      os << "Maximum value = " << *std::max_element(ab.axis_.begin(), ab.axis_.begin() + ab.len_[i]) << std::endl;
    } else {
      os << "Minimum value = " << *std::min_element(ab.axis_.begin() + ab.len_[i - 1], 
              ab.axis_.begin() + ab.len_[i]) << std::endl;
      os << "Maximum value = " << *std::max_element(ab.axis_.begin() + ab.len_[i - 1], 
              ab.axis_.begin() + ab.len_[i]) << std::endl;
    }
  }
  os << "No. of kcoeff = " << ab.kcoeff_.size() << std::endl;
  os << "Minimum value = " << *std::min_element(ab.kcoeff_.begin(), ab.kcoeff_.end()) << std::endl;
  os << "Maximum value = " << *std::max_element(ab.kcoeff_.begin(), ab.kcoeff_.end()) << std::endl;
  return os;
}

Real HitranAbsorber::RefTemp_(Real pres) const {
  int nlevel = refatm_.GetDim1();
  int jl = -1, ju = nlevel, jm;
  // make sure pressure is in ascending order
  while (ju - jl > 1) {
    jm = (ju + jl) >> 1;
    if (pres < refatm_(IPR,jm)) ju = jm;
    else jl = jm;
  }

  // prevent interpolation problem at boundary
  if (jl == -1) jl = 0;
  if (ju == nlevel) ju = nlevel - 1;
  if (jl == ju) return refatm_(IDN,jl);

  //assert(jl >= 0 && ju < nlevel);
  Real result = log(refatm_(IPR,jl)/pres)*log(refatm_(IDN,ju)) 
    + log(pres/refatm_(IPR,ju))*log(refatm_(IDN,jl));
  result = exp(result/log(refatm_(IPR,jl)/refatm_(IPR,ju)));
  return result;
}

void HitranAbsorber::loadCoefficient(std::string fname, int bid)
{
#ifdef NETCDFOUTPUT
  int fileid, dimid, varid, err;
  char tname[80];
  nc_open(fname.c_str(), NC_NETCDF4, &fileid);

  nc_inq_dimid(fileid, "Wavenumber", &dimid);
  nc_inq_dimlen(fileid, dimid, (size_t*)len_);
  nc_inq_dimid(fileid, "Pressure", &dimid);
  nc_inq_dimlen(fileid, dimid, (size_t*)(len_ + 1));
  strcpy(tname, "T_");
  strcat(tname, myname.c_str());
  nc_inq_dimid(fileid, tname, &dimid);
  nc_inq_dimlen(fileid, dimid, (size_t*)(len_ + 2));

  axis_.resize(len_[0] + len_[1] + len_[2]);

  nc_inq_varid(fileid, "Wavenumber", &varid);
  nc_get_var_double(fileid, varid, axis_.data());
  err = nc_inq_varid(fileid, "Pressure", &varid);
  if (err != NC_NOERR)
    throw std::runtime_error(nc_strerror(err));
  err = nc_get_var_double(fileid, varid, axis_.data() + len_[0]);
  if (err != NC_NOERR)
    throw std::runtime_error(nc_strerror(err));
  nc_inq_varid(fileid, tname, &varid);
  nc_get_var_double(fileid, varid, axis_.data() + len_[0] + len_[1]);

  Real *temp = new Real[len_[1]];
  nc_inq_varid(fileid, "Temperature", &varid);
  nc_get_var_double(fileid, varid, temp);

  refatm_.NewAthenaArray(NHYDRO, len_[1]);
  for (int i = 0; i < len_[1]; i++) {
    refatm_(IPR,i) = axis_[len_[0] + i];
    refatm_(IDN,i)  = temp[i];
  }

  kcoeff_.resize(len_[0]*len_[1]*len_[2]);
  nc_inq_varid(fileid, myname.c_str(), &varid);
  nc_get_var_double(fileid, varid, kcoeff_.data());
  nc_close(fileid);
  delete[] temp;
#endif
}

Real HitranAbsorber::getAttenuation(Real wave1, Real wave2,
    CellVariables const& var) const
{
  // first axis is wavenumber, second is pressure, third is temperature anomaly
  Real val, coord[3] = {wave1, var.q[IPR], var.q[IDN] - RefTemp_(var.q[IPR])};
  interpn(&val, coord, kcoeff_.data(), axis_.data(), len_, 3);

  static const Real Rgas = 8.314462;
  Real dens = var.q[IPR]/(Rgas*var.q[IDN]);
  Real x0 = 1.;
  if (imol_ == 0) 
    for (int n = 1; n <= NVAPOR; ++n) x0 -= var.q[n];
  else
    x0 = var.q[imol_];
  return 1.E-3*exp(val)*dens*x0*mixr_; // ln(m*2/kmol) -> 1/m
}
