#pragma once

#include <cmath>
#include <limits>
#include <vector>
#include <algorithm>
#include <stdexcept>

namespace ant {
namespace std_ext {

constexpr double inf = std::numeric_limits<double>::infinity();
constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

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

template<typename T>
constexpr inline T abs_diff(T a, T b) noexcept {
  return a > b ? a - b : b - a;
}

template<class T>
struct RMS_t {
    unsigned n = 0;
    T sum{};
    T sum2{};
    void Add(const T& v) {
        ++n;
        sum += v;
        sum2 += std_ext::sqr(v);
    }
    T GetMean() const {
        return sum/n;
    }
    T GetRMS() const {
        return std::sqrt(std::abs(sum2/n - std_ext::sqr(GetMean())));
    }
    T GetSigmaMean() const {
        return GetRMS() / sqrt(n);
    }
};

using RMS = RMS_t<double>;


struct IQR {
    mutable std::vector<double> nums;
    mutable bool is_sorted = false;

    void Add(double v) {
        nums.emplace_back(v);
        is_sorted = false;
    }

    size_t GetN() const {
        return nums.size();
    }

    double GetMedian() const {
        return GetPercentile(0.5);
    }

    double GetIQR() const {
        return GetPercentile(0.75)-GetPercentile(0.25);
    }


    double GetIQRStdDev() const {
        // scale IQR to match StdDev=sigma=RMS of normal distribution
        // see https://en.wikipedia.org/wiki/Interquartile_range
        return GetIQR()/1.349;
    }

private:
    double GetPercentile(double p) const {
        EnsureSorted();

        if(nums.empty())
            throw std::out_of_range("IQR needs at least one entry");
        const double  x = p*nums.size();

        // approach the right index from the left
        // to handle median correctly
        const double x_l = x - 0.5;
        const double x_r = x + 0.5;
        const unsigned i_l = std::floor(x_l);
        const unsigned i_r = x_r > std::floor(x_r) ? std::floor(x_r) : std::floor(x_r) - 1;

        // calculate the fractional percentile
        return 0.5*nums.at(i_l) + 0.5*nums.at(i_r);
    }

    void EnsureSorted() const {
        if(is_sorted)
            return;
        std::sort(nums.begin(), nums.end());
        is_sorted = true;
    }
};



}} // namespace ant::std_ext

