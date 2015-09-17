#pragma once

#include <cmath>
#include <limits>

namespace ant {
namespace std_ext {

constexpr double inf = std::numeric_limits<double>::infinity();

constexpr inline double degree_to_radian(const double degree) noexcept {
  return degree * M_PI / 180.0;
}
constexpr inline double radian_to_degree(const double radian) noexcept {
  return radian * 180.0 / M_PI;
}

template <typename T>
constexpr inline T sqr(const T& x) noexcept { return x*x; }

}} // namespace ant::std_ext

