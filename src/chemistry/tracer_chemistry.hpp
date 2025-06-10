/** @file tracer_chemistry.hpp
 * @brief
 *
 * @author Ali Hyder (ali.hyder@jpl.nasa.gov)
 * @date Monday 09/06/2025
 * @bug Still testing!!!
 */


#ifndef TRACER_CHEMISTRY_HPP
#define TRACER_CHEMISTRY_HPP

#include <vector>
#include <algorithm>
#include <cmath>


namespace tracerchem {

// Basic chemical functionality to generate thermodynamic properties:
// Using NASA-9 Polynomials...
double EnthalpyPerRT(const std::vector<double> &a, double T);
double EntropyPerR(const std::vector<double> &a, double T);
double GibbsPerRT(const std::vector<double> &a, double T);
// Using NASA-7 Polynomials...
double EnthalpyPerRT_NASA7(const std::vector<double> &a, double T);
double EntropyPerR_NASA7(const std::vector<double> &a, double T);
double GibbsPerRT_NASA7(const std::vector<double> &a, double T);

// Equilibrium rate constants:
double KeqCO(double T);
double KeqPH3(double T);
double KeqPO2(double T);
double KeqGeH2(double T);

// Equilibrium function headers:
double X_CO_Eq(double press, double Temp, double X_CH4, double X_H2O, double X_H2);
double X_PH3_Eq(double Temp, double X_P, double X_H2O, double X_H2);

// Interpolation function for when equilibria are estimated using published sources.
double linearInterpolate(const std::vector<double>& X, const std::vector<double>& T, double x);

}

#endif // TRACER_CHEMISTRY_HPP