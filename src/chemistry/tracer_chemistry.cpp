/** @file tracer_chemistry.cpp
 * @brief
 *
 * @author Ali Hyder (ali.hyder@jpl.nasa.gov)
 * @date Monday 09/06/2025
 * @bug Still testing!!!
 */

#include "tracer_chemistry.hpp"
#include "tracer_chem_NASA_poly.hpp"

#include <cmath>

using namespace NASA_POLY;

namespace tracerchem {

double EnthalpyPerRT(const std::vector<double> &a, double T){
  return -a[0]/pow(T,2) + a[1]*log(T)/T + a[2] + a[3]*T/2 + (a[4]*pow(T,2))/3 + (a[5]*pow(T,3))/4 + (a[6]*pow(T,4))/5 + a[7]/T;
}

double EntropyPerR(const std::vector<double> &a, double T){
  return -(a[0]/pow(T,2))/2 - a[1]/T + a[2]*log(T) + a[3]*T + (a[4]*pow(T,2))/2 + (a[5]*pow(T,3))/3 + (a[6]*pow(T,4))/4 + a[8];                        
}

double EnthalpyPerRT_NASA7(const std::vector<double> &a, double T){
  return a[0] + a[1]*T/2 + (a[2]*pow(T,2))/3 + (a[3]*pow(T,3))/4 + (a[4]*pow(T,4))/5 + a[5]/T;
}

double EntropyPerR_NASA7(const std::vector<double> &a, double T){
  return a[0]*log(T) + a[1]*T + (a[2]*pow(T,2))/2 + (a[3]*pow(T,3))/3 + (a[4]*pow(T,4))/4 + a[6];
}

double GibbsPerRT(const std::vector<double> &a, double T){
  return EnthalpyPerRT(a, T) - EntropyPerR(a, T); // unitless
}

double GibbsPerRT_NASA7(const std::vector<double> &a, double T){
  return EnthalpyPerRT_NASA7(a, T) - EntropyPerR_NASA7(a, T); // unitless
}

double KeqCO(double T){
  
  double delta_Gibbs_f_CO  = GibbsPerRT(GetCoefficients(Species::CO, T), T);
  double delta_Gibbs_f_CH4 = GibbsPerRT(GetCoefficients(Species::CH4, T), T);
  double delta_Gibbs_f_H2O = GibbsPerRT(GetCoefficients(Species::H2O, T), T);

  return exp(-(delta_Gibbs_f_CO - delta_Gibbs_f_CH4 - delta_Gibbs_f_H2O));
}

double KeqPH3(double T){

  double Gibbs_f_H3PO4  = GibbsPerRT_NASA7(GetCoefficients(Species::H3PO4, T), T);
  double Gibbs_f_H2     = GibbsPerRT_NASA7(GetCoefficients(Species::H2, T), T);
  double Gibbs_f_PH3    = GibbsPerRT_NASA7(GetCoefficients(Species::PH3, T), T);
  double Gibbs_f_H2O    = GibbsPerRT_NASA7(GetCoefficients(Species::H2O_N7, T), T);

  return exp(-(Gibbs_f_H3PO4 + 4*Gibbs_f_H2 - Gibbs_f_PH3 - 4*Gibbs_f_H2O));
}

double KeqPO2(double T){

  double Gibbs_f_PO2   = GibbsPerRT_NASA7(GetCoefficients(Species::PO2, T), T);
  double Gibbs_f_H     = GibbsPerRT_NASA7(GetCoefficients(Species::H, T), T);
  double Gibbs_f_HOPO2 = GibbsPerRT_NASA7(GetCoefficients(Species::HOPO2, T), T);
  double Gibbs_f_H2O   = GibbsPerRT_NASA7(GetCoefficients(Species::H2O_N7, T), T);

  return exp(-(Gibbs_f_HOPO2 + Gibbs_f_H - Gibbs_f_PO2 - Gibbs_f_H2O));
}

double KeqGeH2(double T){
  
  double Gibbs_f_GeH2 = GibbsPerRT_NASA7(GetCoefficients(Species::GeH2, T), T);
  double Gibbs_f_H2   = GibbsPerRT_NASA7(GetCoefficients(Species::H2, T), T);
  double Gibbs_f_GeH4 = GibbsPerRT_NASA7(GetCoefficients(Species::GeH4, T), T);

  return exp(-(Gibbs_f_GeH2 + Gibbs_f_H2 - Gibbs_f_GeH4));
}

double X_CO_Eq(double press, double Temp, double X_CH4, double X_H2O, double X_H2){

  double K_eq_CO = KeqCO(Temp);
  return X_CH4 * X_H2O * K_eq_CO / pow(X_H2, 3) / pow((press/1E5), 2);

}

double X_PH3_Eq(double Temp, double X_P, double X_H2O, double X_H2){

  double K_eq_PH3 = KeqPH3(Temp);
  return X_P / (1 + (K_eq_PH3 * pow(X_H2O/X_H2, 4)));

}

// Function to perform linear interpolation.
// Given a value x and arrays X and T, it finds the appropriate interval in X and interpolates the corresponding value in T.
double linearInterpolate(const std::vector<double>& X, const std::vector<double>& T, double x) {
    // Check if x is out of bounds of the X array
    if (x <= X.front()) return T.front();
    if (x >= X.back()) return T.back();
    // Find the right place in the array (lower bound of the interval)
    auto lower = std::lower_bound(X.begin(), X.end(), x);
    int index = lower - X.begin();
    // Edge case handling: if x matches the last element, interpolate with the second last
    if (index == X.size()) {
        index--;
    }
    // Calculate the interpolation
    double x1 = X[index - 1], x2 = X[index];
    double t1 = T[index - 1], t2 = T[index];
    // The interpolation formula
    double t = t1 + (x - x1) * (t2 - t1) / (x2 - x1);
    return t;
}

}