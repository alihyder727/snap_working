// Athena++ header
#include "../thermodynamics/thermodynamics.hpp"
#include "hydro.hpp"

void Hydro::VaporConcentrationLimiter(AthenaArray<Real> &u)
{
  MeshBlock *pmb = pmy_block;
  Thermodynamics *pthermo = pmb->pthermo;
  if (NVAPOR == 0) return;

  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  for (int k=ks; k<=ke; ++k)
    for (int j=js; j<=je; ++j) {
      // fix negative vapor from top down
      for (int i = ie; i > is; --i) {
        Real density = 0., v1, v2, v3;

        for (int n = 0; n <= NVAPOR; ++n) density += u(n,k,j,i);
        v1 = u(IM1,k,j,i)/density,
        v2 = u(IM2,k,j,i)/density,
        v3 = u(IM3,k,j,i)/density;

        // internal energy
        Real KE = 0.5*density*(v1*v1 + v2*v2 + v3*v3);
        Real LE = 0., rhocv = 0.;
        for (int n = 0; n <= NVAPOR; ++n)
          rhocv += u(n,k,j,i)*pthermo->GetCv(n);
        Real temp = (u(IEN,k,j,i) - KE)/rhocv;

        for (int n = 1; n <= NVAPOR; ++n) {
          Real rho = u(n,k,j,i);
          if (rho < 0.) {
            u(n,k,j,i-1) += rho;
            u(IM1,k,j,i-1) += v1*rho;
            u(IM2,k,j,i-1) += v2*rho;
            u(IM3,k,j,i-1) += v3*rho;
            Real en = pthermo->GetCv(n)*temp + 0.5*(v1*v1 + v2*v2 + v3*v3);
            u(IEN,k,j,i-1) += en*rho;

            u(n,k,j,i) = 0.;
            u(IM1,k,j,i) -= v1*rho;
            u(IM2,k,j,i) -= v2*rho;
            u(IM3,k,j,i) -= v3*rho;
            u(IEN,k,j,i) -= en*rho;
          }
        }
      }
    }
}
