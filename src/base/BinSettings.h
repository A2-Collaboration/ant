#pragma once

#include "base/interval.h"
#include <vector>
#include <stdexcept>

namespace ant {

class BinSettings: public interval<double>  {
protected:
    unsigned int bins;
public:

    struct Exception : std::runtime_error {
        using std::runtime_error::runtime_error;
    };

    BinSettings(unsigned int number_of_bins, double minimum, double maximum):
        interval<double>(minimum,maximum),
        bins(number_of_bins)
    { if(!IsSane()) throw Exception("Max < Min"); }

    BinSettings(unsigned int number_of_bins):
        interval<double>(0,number_of_bins),
        bins(number_of_bins)
    { if(!IsSane()) throw Exception("Max < Min"); }

    BinSettings(unsigned int number_of_bins, const interval<double>& i):
        interval<double>(i),
        bins(number_of_bins)
    { if(!IsSane()) throw Exception("Max < Min"); }

    unsigned int Bins() const noexcept { return bins; }
    unsigned int& Bins()   noexcept { return bins; }

    double BinWidth() const { return Length() / bins; }

    /**
     * @brief Produces BinSettings close to the supplied one, but with a given bin size.
     *        Slightly modifies the number of bins and boundaries.
     *        Useful to to counteract binning effects if intrinsic spacing of data is known (for example ADC channel widths).
     * @param bins Target settings
     * @param binSize size to match
     * @return new settings with the given bin size
     */
    static BinSettings RoundToBinSize(const BinSettings& bins, const double binSize);

    /**
     * @brief Make creates optimal BinSettings for given x_values
     * @param x_values strictly monotonic list of possible values
     * @return optimal BinSettings
     */
    static BinSettings Make(const std::vector<double>& x_values);

    /**
     * @brief Get the number of the bin a value would end up in
     * @param v The Value
     * @return a bin number between 0 (lowest bin) to BinSettings::Bins-1, -1 if v is outside the interval
     * @note Bin numbering is different from ROOTs TH{1..3} numbering, where the lowest bin has the index 1
     */
    int getBin(const double v) const noexcept;
};

}
