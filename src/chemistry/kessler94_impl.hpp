/** @file kessler94_impl.hpp
 * @brief Implement Kessler94 Chemistry
 *
 * This file is automatically generated by make_chemistry.py
 *
 * @author Cheng Li
 * @bug No know bugs.
 */

template<typename D1, typename D2>
void Kessler94::AssembleReactionMatrix(Eigen::DenseBase<D1>& rate,
  Eigen::DenseBase<D2>& jac, Real const q[], Real cv, Real time)
{
  Thermodynamics *pthermo = pmy_chem->pmy_block->pthermo;

  int iT = index_[0];
  Real T = q[iT];
  int iqv = index_[1];
  Real qv = q[iqv];
  int iqc = index_[2];
  Real qc = q[iqc];
  int iqp = index_[3];
  Real qp = q[iqp];

  pthermo->SaturationSurplus(dqsat_.data(), q, VariableType::chem);
  Real dq = dqsat_[iqv];
  Real qs = qv - dq;
  Real Rv = pthermo->GetRd()/pthermo->GetMassRatio(iqv);
  Real lf = pthermo->GetLatent(NVAPOR+iqv,T) - Rv*T;
  Real dqsdt = qs/T*lf/(Rv*T);
  Real mols = q[IPR]/(Thermodynamics::Rgas*q[IDN]);

  Real k1 = coeffs_["condensation"];
  Real k2 = coeffs_["autoconversion"];
  Real k3 = coeffs_["accretion"];
  Real k4 = coeffs_["evaporation"];

  // qc -> qp; k2
  rate(2) -= k2*qc;
  jac(2,2) -= k2;
  
  rate(3) += k2*qc;
  jac(3,2) += k2;
  
  // enthalpy change
  rate(0) += (k2*qc)*(deltaU_[iqc] - deltaU_[iqp])/cv;
  jac(0,2) += (k2)*(deltaU_[iqc] - deltaU_[iqp])/cv;
  
  // qc + qp -> 2qp; k3
  rate(2) -= k3*qc*qp;
  jac(2,2) -= k3*qp;
  jac(2,3) -= k3*qc;
  
  rate(3) += k3*qc*qp;
  jac(3,2) += k3*qp;
  jac(3,3) += k3*qc;
  
  // enthalpy change
  rate(0) += (k3*qc*qp)*(deltaU_[iqc] - deltaU_[iqp])/cv;
  jac(0,2) += (k3*qp)*(deltaU_[iqc] - deltaU_[iqp])/cv;
  jac(0,3) += (k3*qc)*(deltaU_[iqc] - deltaU_[iqp])/cv;
  
  if (qv < qs) {
    // qc -> qv; k1*(qs - qv)/mols
    rate(2) -= k1*(qs - qv)/mols*qc;
    jac(2,1) -= -k1*qc/mols;
    jac(2,2) -= k1*(qs - qv)/mols;
    
    rate(1) += k1*(qs - qv)/mols*qc;
    jac(1,1) += -k1*qc/mols;
    jac(1,2) += k1*(qs - qv)/mols;
    
    // enthalpy change
    rate(0) += (k1*(qs - qv)/mols*qc)*(deltaU_[iqc] - deltaU_[iqv])/cv;
    jac(0,1) += (-k1*qc/mols)*(deltaU_[iqc] - deltaU_[iqv])/cv;
    jac(0,2) += (k1*(qs - qv)/mols)*(deltaU_[iqc] - deltaU_[iqv])/cv;
    
    // qp -> qv; k4*(qs - qv)/mols
    rate(3) -= k4*(qs - qv)/mols*qp;
    jac(3,1) -= -k4*qp/mols;
    jac(3,3) -= k4*(qs - qv)/mols;
    
    rate(1) += k4*(qs - qv)/mols*qp;
    jac(1,1) += -k4*qp/mols;
    jac(1,3) += k4*(qs - qv)/mols;
    
    // enthalpy change
    rate(0) += (k4*(qs - qv)/mols*qp)*(deltaU_[iqp] - deltaU_[iqv])/cv;
    jac(0,1) += (-k4*qp/mols)*(deltaU_[iqp] - deltaU_[iqv])/cv;
    jac(0,3) += (k4*(qs - qv)/mols)*(deltaU_[iqp] - deltaU_[iqv])/cv;
  }
  
  if (qv > qs) {
    // qv -> qc; k1*(qv - qs)
    rate(1) -= k1*(qv - qs);
    jac(1,0) -= -k1*dqsdt;
    jac(1,1) -= k1;
    
    rate(2) += k1*(qv - qs);
    jac(2,0) += -k1*dqsdt;
    jac(2,1) += k1;
    
    // enthalpy change
    rate(0) += (k1*(qv - qs))*(deltaU_[iqv] - deltaU_[iqc])/cv;
    jac(0,0) += (-k1*dqsdt)*(deltaU_[iqv] - deltaU_[iqc])/cv;
    jac(0,1) += (k1)*(deltaU_[iqv] - deltaU_[iqc])/cv;
  }

}
