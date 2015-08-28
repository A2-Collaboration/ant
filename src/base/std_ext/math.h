#pragma once

#include <cmath>
#include <limits>

namespace ant {
namespace std_ext {

constexpr double inf = std::numeric_limits<double>::infinity();

inline double degree_to_radian(const double degree) {
  return degree * M_PI / 180.0;
}
inline double radian_to_degree(const double radian) {
  return radian * 180.0 / M_PI;
}

template <typename T>
inline T sqr(const T& x) { return x*x; }

}} // namespace ant::std_ext

