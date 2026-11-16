#ifndef HARVEST_MATH_HPP
#define HARVEST_MATH_HPP
#include <cmath>
#include <stdexcept>
#include <limits>
#include <vector>
namespace harvest_h03 {
struct Rule { double lo, hi, fraction; };
inline double removal_fraction(double age, const std::vector<Rule>& rules) {
  for (const auto& rule : rules) {
    if (age >= rule.lo && age <= rule.hi) return rule.fraction;
  }
  return 0.0;
}

inline double removed(double biomass, double fraction) {
  if (!std::isfinite(biomass) || biomass < 0 || biomass != std::trunc(biomass) || biomass > 2147483647.0)
    throw std::invalid_argument("H03-A: biomass must be an integer in [0, INT32_MAX]");
  if (!std::isfinite(fraction) || fraction < 0 || fraction > 1)
    throw std::invalid_argument("H03-A: fraction must be in [0,1]");
  return std::trunc(biomass * fraction);
}
}
#endif
