#ifndef FREEDMAN_SIMPLE_HPP
#define FREEDMAN_SIMPLE_HPP

// C/C++ header
#include <vector>

// Athena++ header
#include "absorber.hpp"

class FreedmanSimple: public Absorber {
public:
  FreedmanSimple(RadiationBand *pband, ParameterInput *pin);
  virtual ~FreedmanSimple() {}
  Real getAttenuation(Real wave1, Real wave2, CellVariables const& var) const;

private:
  Real scale_;
};

#endif
