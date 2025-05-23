//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file globals.cpp
//  \brief namespace containing global variables.
//
// Yes, we all know global variables should NEVER be used, but in fact they are ideal for,
// e.g., global constants that are set once and never changed.  To prevent name collisions
// global variables are wrapped in their own namespace.

// C headers

// C++ headers

// Athena++ headers
#include "athena.hpp"
#include "globals.hpp"

namespace Globals {
// all of these global variables are set at the start of main():
int my_rank;         // MPI rank of this process
int nranks;          // total number of MPI ranks
int mpi_tag_ub;
}
