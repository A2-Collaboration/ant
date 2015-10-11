#pragma once

#include <cmath>
#include <limits>

namespace ant {
namespace std_ext {

constexpr double inf = std::numeric_limits<double>::infinity();

template <typename T>
constexpr inline T degree_to_radian(const T& degree) noexcept {
  return degree * M_PI / 180.0;
}

template <typename T>
constexpr inline T radian_to_degree(const T& radian) noexcept {
  return radian * 180.0 / M_PI;
}

template <typename T>
constexpr inline T sqr(const T& x) noexcept { return x*x; }

}} // namespace ant::std_ext

