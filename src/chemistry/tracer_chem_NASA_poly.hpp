/** @file tracer_chem_NASA_poly.hpp
 * @brief
 *
 * @author Ali Hyder (ali.hyder@jpl.nasa.gov)
 * @date Monday 09/06/2025
 * @bug Still testing!!!
 */

#ifndef TRACER_CHEM_NASA_POLY_HPP
#define TRACER_CHEM_NASA_POLY_HPP

#include <vector>
#include <string>
#include <unordered_map>

namespace NASA_POLY {

enum class Species {
    CO    = 0,
    CH4   = 1,
    H2O   = 2,
    H3PO4 = 3,
    H2O_N7= 4,
    H2    = 5,
    H     = 6,
    PO2   = 7,
    HOPO2 = 8,
    GeH4  = 9,
    GeH2  = 10,
    PH3   = 11
};

// Structure to hold low/high coefficients and switching temperature
struct NASACoeffPair {
  std::vector<double> low;
  std::vector<double> high;
  double T_mid;  // Typically 1000 K in NASA data
};

// Main species data map
// extern const std::unordered_map<std::string, NASACoeffPair> species_data;

// Function that chooses the right coefficient set based on temperature

// Using unordered map - similar to enum but it may be more time consuming...
// const std::vector<double>& GetCoefficients(const std::string &species, double T);

// Using enum-based species selection:
const std::vector<double>& GetCoefficients(Species sp, double T);

} // namespace NASA_POLY

#endif // TRACER_CHEM_NASA_POLY_HPP
