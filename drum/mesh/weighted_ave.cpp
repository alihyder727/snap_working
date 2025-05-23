//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file weighted_ave.cpp
//  \brief

// C headers

// C++ headers

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../debugger/debugger.hpp"
#include "mesh.hpp"

//----------------------------------------------------------------------------------------
//! \fn  void WeightedAve::WeightedAve
//  \brief Compute weighted average of AthenaArrays (including cell-averaged U in time
//         integrator step)

void MeshBlock::WeightedAve(AthenaArray<Real> &u_out, AthenaArray<Real> &u_in1,
                            AthenaArray<Real> &u_in2, const Real wght[3]) {
  // consider every possible simplified form of weighted sum operator:
  // U = a*U + b*U1 + c*U2

  // assuming all 3x arrays are of the same size (or at least u_out is equal or larger
  // than each input array) in each array dimension, and full range is desired:
  // nx4*(3D real MeshBlock cells)
  const int nu = u_out.GetDim4() - 1;
  pdebug->Call("MeshBlock::WeightedAvg");

  // u_in2 may be an unallocated AthenaArray if using a 2S time integrator
  if (wght[0] == 1.0) {
    if (wght[2] != 0.0) {
      for (int n=0; n<=nu; ++n) {
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
#pragma omp simd
            for (int i=is; i<=ie; ++i) {
              u_out(n,k,j,i) += wght[1]*u_in1(n,k,j,i) + wght[2]*u_in2(n,k,j,i);
            }
          }
        }
      }
    } else { // do not dereference u_in2
      if (wght[1] != 0.0) {
        for (int n=0; n<=nu; ++n) {
          for (int k=ks; k<=ke; ++k) {
            for (int j=js; j<=je; ++j) {
#pragma omp simd
              for (int i=is; i<=ie; ++i) {
                u_out(n,k,j,i) += wght[1]*u_in1(n,k,j,i);
              }
            }
          }
        }
      }
    }
  } else if (wght[0] == 0.0) {
    if (wght[2] != 0.0) {
      for (int n=0; n<=nu; ++n) {
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
#pragma omp simd
            for (int i=is; i<=ie; ++i) {
              u_out(n,k,j,i) = wght[1]*u_in1(n,k,j,i) + wght[2]*u_in2(n,k,j,i);
            }
          }
        }
      }
    } else if (wght[1] == 1.0) {
      // just deep copy
      for (int n=0; n<=nu; ++n) {
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
#pragma omp simd
            for (int i=is; i<=ie; ++i) {
              u_out(n,k,j,i) = u_in1(n,k,j,i);
            }
          }
        }
      }
    } else {
      for (int n=0; n<=nu; ++n) {
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
#pragma omp simd
            for (int i=is; i<=ie; ++i) {
              u_out(n,k,j,i) = wght[1]*u_in1(n,k,j,i);
            }
          }
        }
      }
    }
  } else {
    if (wght[2] != 0.0) {
      for (int n=0; n<=nu; ++n) {
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
#pragma omp simd
            for (int i=is; i<=ie; ++i) {
              u_out(n,k,j,i) = wght[0]*u_out(n,k,j,i) + wght[1]*u_in1(n,k,j,i)
                               + wght[2]*u_in2(n,k,j,i);
            }
          }
        }
      }
    } else { // do not dereference u_in2
      if (wght[1] != 0.0) {
        for (int n=0; n<=nu; ++n) {
          for (int k=ks; k<=ke; ++k) {
            for (int j=js; j<=je; ++j) {
#pragma omp simd
              for (int i=is; i<=ie; ++i) {
                u_out(n,k,j,i) = wght[0]*u_out(n,k,j,i) + wght[1]*u_in1(n,k,j,i);
              }
            }
          }
        }
      } else { // do not dereference u_in1
        for (int n=0; n<=nu; ++n) {
          for (int k=ks; k<=ke; ++k) {
            for (int j=js; j<=je; ++j) {
#pragma omp simd
              for (int i=is; i<=ie; ++i) {
                u_out(n,k,j,i) *= wght[0];
              }
            }
          }
        }
      }
    }
  }

#if DEBUG_LEVEL > 2
  pdebug->CheckConservation("u_out", u_out, is, ie, js, je, ks, ke);
#endif

  pdebug->Leave();
  return;
}


//----------------------------------------------------------------------------------------
//! \fn  void MeshBlock::WeightedAve
//  \brief Compute weighted average of face-averaged B in time integrator step

void MeshBlock::WeightedAve(FaceField &b_out, FaceField &b_in1, FaceField &b_in2,
                            const Real wght[3]) {
  int jl=js; int ju=je+1;
  // move these limit modifications outside the loop
  if (pbval->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar
      || pbval->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar_wedge)
    jl=js+1;
  if (pbval->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar
      || pbval->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar_wedge)
    ju=je;

  // Note: these loops can be combined now that they avoid curl terms
  // Only need to separately account for the final longitudinal face in each loop limit
  if (wght[0] == 1.0) {
    if (wght[2] != 0.0) {
      //---- B1
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma omp simd
          for (int i=is; i<=ie+1; ++i) {
            b_out.x1f(k,j,i) += wght[1]*b_in1.x1f(k,j,i) + wght[2]*b_in2.x1f(k,j,i);
          }
        }
      }
      //---- B2
      for (int k=ks; k<=ke; ++k) {
        for (int j=jl; j<=ju; ++j) {
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            b_out.x2f(k,j,i) += wght[1]*b_in1.x2f(k,j,i) + wght[2]*b_in2.x2f(k,j,i);
          }
        }
      }
      //---- B3
      for (int k=ks; k<=ke+1; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            b_out.x3f(k,j,i) += wght[1]*b_in1.x3f(k,j,i) + wght[2]*b_in2.x3f(k,j,i);
          }
        }
      }
    } else { // do not dereference u_in2
      if (wght[1] != 0.0) {
        //---- B1
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
#pragma omp simd
            for (int i=is; i<=ie+1; ++i) {
              b_out.x1f(k,j,i) += wght[1]*b_in1.x1f(k,j,i);
            }
          }
        }
        //---- B2
        for (int k=ks; k<=ke; ++k) {
          for (int j=jl; j<=ju; ++j) {
#pragma omp simd
            for (int i=is; i<=ie; ++i) {
              b_out.x2f(k,j,i) += wght[1]*b_in1.x2f(k,j,i);
            }
          }
        }
        //---- B3
        for (int k=ks; k<=ke+1; ++k) {
          for (int j=js; j<=je; ++j) {
#pragma omp simd
            for (int i=is; i<=ie; ++i) {
              b_out.x3f(k,j,i) += wght[1]*b_in1.x3f(k,j,i);
            }
          }
        }
      }
    }
  } else if (wght[0] == 0.0) {
    if (wght[2] != 0.0) {
      //---- B1
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma omp simd
          for (int i=is; i<=ie+1; ++i) {
            b_out.x1f(k,j,i) = wght[1]*b_in1.x1f(k,j,i) + wght[2]*b_in2.x1f(k,j,i);
          }
        }
      }
      //---- B2
      for (int k=ks; k<=ke; ++k) {
        for (int j=jl; j<=ju; ++j) {
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            b_out.x2f(k,j,i) = wght[1]*b_in1.x2f(k,j,i) + wght[2]*b_in2.x2f(k,j,i);
          }
        }
      }
      //---- B3
      for (int k=ks; k<=ke+1; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            b_out.x3f(k,j,i) = wght[1]*b_in1.x3f(k,j,i) + wght[2]*b_in2.x3f(k,j,i);
          }
        }
      }
    } else if (wght[1] == 1.0) {
      // just deep copy
      //---- B1
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma omp simd
          for (int i=is; i<=ie+1; ++i) {
            b_out.x1f(k,j,i) = b_in1.x1f(k,j,i);
          }
        }
      }
      //---- B2
      for (int k=ks; k<=ke; ++k) {
        for (int j=jl; j<=ju; ++j) {
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            b_out.x2f(k,j,i) = b_in1.x2f(k,j,i);
          }
        }
      }
      //---- B3
      for (int k=ks; k<=ke+1; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            b_out.x3f(k,j,i) = b_in1.x3f(k,j,i);
          }
        }
      }
    } else {
      //---- B1
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma omp simd
          for (int i=is; i<=ie+1; ++i) {
            b_out.x1f(k,j,i) = wght[1]*b_in1.x1f(k,j,i);
          }
        }
      }
      //---- B2
      for (int k=ks; k<=ke; ++k) {
        for (int j=jl; j<=ju; ++j) {
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            b_out.x2f(k,j,i) = wght[1]*b_in1.x2f(k,j,i);
          }
        }
      }
      //---- B3
      for (int k=ks; k<=ke+1; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            b_out.x3f(k,j,i) = wght[1]*b_in1.x3f(k,j,i);
          }
        }
      }
    }
  } else {
    if (wght[2] != 0.0) {
      //---- B1
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma omp simd
          for (int i=is; i<=ie+1; ++i) {
            b_out.x1f(k,j,i) = wght[0]*b_out.x1f(k,j,i) + wght[1]*b_in1.x1f(k,j,i)
                               + wght[2]*b_in2.x1f(k,j,i);
          }
        }
      }
      //---- B2
      for (int k=ks; k<=ke; ++k) {
        for (int j=jl; j<=ju; ++j) {
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            b_out.x2f(k,j,i) = wght[0]*b_out.x2f(k,j,i) + wght[1]*b_in1.x2f(k,j,i)
                               + wght[2]*b_in2.x2f(k,j,i);
          }
        }
      }
      //---- B3
      for (int k=ks; k<=ke+1; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            b_out.x3f(k,j,i) = wght[0]*b_out.x3f(k,j,i) + wght[1]*b_in1.x3f(k,j,i)
                               + wght[2]*b_in2.x3f(k,j,i);
          }
        }
      }
    } else { // do not dereference u_in2
      if (wght[1] != 0.0) {
        //---- B1
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
#pragma omp simd
            for (int i=is; i<=ie+1; ++i) {
              b_out.x1f(k,j,i) = wght[0]*b_out.x1f(k,j,i) + wght[1]*b_in1.x1f(k,j,i);
            }
          }
        }
        //---- B2
        for (int k=ks; k<=ke; ++k) {
          for (int j=jl; j<=ju; ++j) {
#pragma omp simd
            for (int i=is; i<=ie; ++i) {
              b_out.x2f(k,j,i) = wght[0]*b_out.x2f(k,j,i) + wght[1]*b_in1.x2f(k,j,i);
            }
          }
        }
        //---- B3
        for (int k=ks; k<=ke+1; ++k) {
          for (int j=js; j<=je; ++j) {
#pragma omp simd
            for (int i=is; i<=ie; ++i) {
              b_out.x3f(k,j,i) = wght[0]*b_out.x3f(k,j,i) + wght[1]*b_in1.x3f(k,j,i);
            }
          }
        }
      } else { // do not dereference u_in1
        //---- B1
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
#pragma omp simd
            for (int i=is; i<=ie+1; ++i) {
              b_out.x1f(k,j,i) *= wght[0];
            }
          }
        }
        //---- B2
        for (int k=ks; k<=ke; ++k) {
          for (int j=jl; j<=ju; ++j) {
#pragma omp simd
            for (int i=is; i<=ie; ++i) {
              b_out.x2f(k,j,i) *= wght[0];
            }
          }
        }
        //---- B3
        for (int k=ks; k<=ke+1; ++k) {
          for (int j=js; j<=je; ++j) {
#pragma omp simd
            for (int i=is; i<=ie; ++i) {
              b_out.x3f(k,j,i) *= wght[0];
            }
          }
        }
      }
    }
  }
  return;
}
