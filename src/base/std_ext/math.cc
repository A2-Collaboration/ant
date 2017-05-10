#include "math.h"

#include <algorithm>
#include <stdexcept>


int ant::std_ext::calcNchooseK(int n, int k) {
    int result = 1;
    for(int i=1;i<=k;i++) {
        result *= n - (k - i);
        result /= i;
    }
    return result;
}

void ant::std_ext::IQR::Add(double v) {
    if(std::isfinite(v))
        nums.emplace_back(v);
    is_sorted = false;
}

std::size_t ant::std_ext::IQR::GetN() const {
    return nums.size();
}

double ant::std_ext::IQR::GetMedian() const {
    return GetPercentile(0.5);
}

double ant::std_ext::IQR::GetIQR() const {
    return GetPercentile(0.75)-GetPercentile(0.25);
}

double ant::std_ext::IQR::GetIQRStdDev() const {
    // scale IQR to match StdDev=sigma=RMS of normal distribution
    // see https://en.wikipedia.org/wiki/Interquartile_range
    return GetIQR()/1.349;
}

double ant::std_ext::IQR::GetPercentile(double p) const {
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

void ant::std_ext::IQR::EnsureSorted() const {
    if(is_sorted)
        return;
    std::sort(nums.begin(), nums.end());
    is_sorted = true;
}
