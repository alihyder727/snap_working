// Athena++ header files
#include "validate_chemistry.hpp"
#include "chemistry_solver.hpp"
#include "../debugger/debugger.hpp"

template<typename T>
void ChemistryBase<T>::IntegrateDense(AthenaArray<Real> &uh, AthenaArray<Real> &up,
  Real time, Real dt)
{
  //std::cout << "Integrating " << myname << std::endl;
  //for (int n = 0; n < 4; ++n)
  //  std::cout << index_[n] << " ";
  //std::cout << std::endl;
  MeshBlock *pmb = pmy_chem->pmy_block;
  pmb->pdebug->Call("Chemistry::IntegrateDense-" + myname);
  std::stringstream &msg = pmb->pdebug->msg;

  Thermodynamics *pthermo = pmb->pthermo;
  EquationOfState *peos = pmb->peos;
  Particles *ppart = pmb->ppart->FindParticle(particle_name);

  int ks = pmb->ks, js = pmb->js, is = pmb->is;
  int ke = pmb->ke, je = pmb->je, ie = pmb->ie;

  int size = uh.GetDim4() + up.GetDim4();

/*  Real total = 0.;
  if (myname == "H2O")
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i)
          total += uh(1,k,j,i) + up(0,k,j,i) + up(1,k,j,i);
  else if (myname == "NH3")
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i)
          total += uh(2,k,j,i) + up(0,k,j,i) + up(1,k,j,i);
#ifdef MPI_PARALLEL
  if (Globals::my_rank == 0) {
    MPI_Reduce(MPI_IN_PLACE, &total, 1, MPI_ATHENA_REAL, MPI_SUM, 0, MPI_COMM_WORLD);
  } else {
    MPI_Reduce(&total, &total, 1, MPI_ATHENA_REAL, MPI_SUM, 0, MPI_COMM_WORLD);
  }
#endif
  msg << "- total before = " << total << std::endl; */

  Real* c0 = new Real [size];
  Real* c1 = new Real [size];
  Real* c2 = new Real [size];

  T* pchem = static_cast<T*>(this);
  Eigen::Matrix<Real, T::Solver::Size, T::Solver::Size> Jac;
  Eigen::Matrix<Real, T::Solver::Size, 1> Rate, sol;

  Real cvd = pthermo->GetRd()/(peos->GetGamma() - 1.);

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        // 1. set rate and jacobian matrix to zero
        Rate.setZero();
        Jac.setZero();

        // 2. from mass density to molar density
        // ignore the kinetic energy of condensates
        Real ux = uh(IVX,k,j,i), uy = uh(IVY,k,j,i), uz = uh(IVZ,k,j,i);
        Real rho_total = 0.;
        for (int n = 1; n <= NVAPOR; ++n)
          rho_total += uh(n,k,j,i);
        for (int n = 0; n < up.GetDim4(); ++n)
          rho_total += up(n,k,j,i);

        pthermo->ConservedToChemical(c0, uh.at(k,j,i));
        Real cd = c0[IPR]/(Thermodynamics::Rgas*c0[IDN]);
        for (int n = 1; n <= NVAPOR; ++n) cd -= c0[n];

        for (int n = 0; n < up.GetDim4(); ++n)
          c0[NHYDRO+n] = up(n,k,j,i)/ppart->GetMolecularWeight(n);

        // 3. make two copies, one for BDF1 and another for TR-BDF2
        std::memcpy(c1, c0, size*sizeof(Real));
        std::memcpy(c2, c0, size*sizeof(Real));

        // 4. calculate heat capacity 
        Real cvt = 0.;
        for (int n = 0; n <= NVAPOR; ++n)
          cvt += uh(n,k,j,i)*pthermo->GetCvRatio(n)*cvd;
        for (int n = 0; n < up.GetDim4(); ++n)
          cvt += up(n,k,j,i)*ppart->GetCv(n);

        // 5. set reaction rate and jacobian matrix
        pchem->AssembleReactionMatrix(Rate, Jac, c0, cvt, time);

        // 6. BDF1 solver
        sol = pchem->solver.BDF1(Rate, Jac, dt);
        for (int n = 0; n < T::Solver::Size; ++n)
          c1[index_[n]] += sol(n);

        /* 7. TR-BDF2 solver
        sol = pchem->solver.TRBDF2(Rate, Jac, dt);
        for (int n = 0; n < T::Solver::Size; ++n)
          c2[index_[n]] += sol(n);

        // 8. Blend solver
        //sol = pchem->solver.TRBDF2_Blend(c2, Rate, Jac, index_.data(), dt);
        Real alpha = 1.;
        for (int n = 1; n < T::Solver::Size; ++n)
          if (c2[index_[n]] < 0.)
            alpha = std::min(alpha, c1[index_[n]]/(c1[index_[n]] - c2[index_[n]]));
        for (int n = 0; n < T::Solver::Size; ++n)
          c2[index_[n]] = (1. - alpha)*c1[index_[n]] + alpha*c2[index_[n]];

        if (c1[1] > c1[IPR]/(c1[0]*Thermodynamics::Rgas)) {
          std::cout << cvt << std::endl;
          for (int n = 0; n < NHYDRO; ++n)
            std::cout << uh(n,k,j,i) << " ";
          for (int n = 0; n < 2; ++n)
            std::cout << up(n,k,j,i) << " ";
          std::cout << std::endl;
          for (int n = 0; n < NHYDRO+2; ++n)
            std::cout << c0[n] << " ";
          std::cout << std::endl;
          for (int n = 0; n < NHYDRO+2; ++n)
            std::cout << c1[n] << " ";
          std::cout << std::endl << std::endl;
          exit(1);
        }*/

        // 9. Apply solution
        //pchem->ApplySolutionLimit(c1, c0, sol, cvt);

        // 10. Adjust pressure
        c1[IPR] = cd*Thermodynamics::Rgas*c1[IDN];
        for (int n = 1; n <= NVAPOR; ++n)
          c1[IPR] += c1[n]*Thermodynamics::Rgas*c1[IDN];

        // 11. from molar density to mass density
        pthermo->ChemicalToConserved(uh.at(k,j,i), c1);
        // reset kinetic energy
        uh(IVX,k,j,i) = ux;
        uh(IVY,k,j,i) = uy;
        uh(IVZ,k,j,i) = uz;

        for (int n = 1; n <= NVAPOR; ++n)
          rho_total -= uh(n,k,j,i);
        Real par_total = 0.;
        for (int n = 0; n < up.GetDim4(); ++n)
          par_total += c1[NHYDRO+n]*ppart->GetMolecularWeight(n);
        if (par_total > 0.)
          for (int n = 0; n < up.GetDim4(); ++n)
            up(n,k,j,i) = c1[NHYDRO+n]*ppart->GetMolecularWeight(n)*rho_total/par_total;
        else
          for (int n = 0; n < up.GetDim4(); ++n)
            up(n,k,j,i) = c1[NHYDRO+n]*ppart->GetMolecularWeight(n);
      }

#if DEBUG_LEVEL > 2
  pmb->pdebug->CheckConservation("uh", uh, is, ie, js, je, ks, ke);
  pmb->pdebug->CheckConservation("up", up, is, ie, js, je, ks, ke);
#endif

/*  total = 0.;
  if (myname == "H2O")
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i)
          total += uh(1,k,j,i) + up(0,k,j,i) + up(1,k,j,i);
  else if (myname == "NH3")
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i)
          total += uh(2,k,j,i) + up(0,k,j,i) + up(1,k,j,i);
#ifdef MPI_PARALLEL
  if (Globals::my_rank == 0) {
    MPI_Reduce(MPI_IN_PLACE, &total, 1, MPI_ATHENA_REAL, MPI_SUM, 0, MPI_COMM_WORLD);
  } else {
    MPI_Reduce(&total, &total, 1, MPI_ATHENA_REAL, MPI_SUM, 0, MPI_COMM_WORLD);
  }
#endif
  msg << "- total after = " << total << std::endl;*/

  delete[] c0;
  delete[] c1;
  delete[] c2;
  pmb->pdebug->Leave();
}
