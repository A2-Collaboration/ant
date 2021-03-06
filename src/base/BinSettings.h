#pragma once

#include "base/interval.h"
#include <vector>
#include <stdexcept>
#include <istream>
#include <ostream>

namespace ant {

class BinSettings : public interval<double>  {
protected:
    unsigned int bins;
public:

    struct Exception : std::runtime_error {
        using std::runtime_error::runtime_error;
    };

    BinSettings(unsigned int number_of_bins, double minimum, double maximum) :
        BinSettings(number_of_bins, {minimum, maximum})
    {}

    BinSettings(unsigned int number_of_bins) :
        BinSettings(number_of_bins, {0, 1.0*number_of_bins})
    {}

    BinSettings(unsigned int number_of_bins, const interval<double>& i):
        interval<double>(i),
        bins(number_of_bins)
    { if(!IsSane()) throw Exception("BinSettings: Max < Min"); }

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

    // the >> operator parses stringified versions made with <<
    friend std::istream& operator>>(std::istream& in, BinSettings& t);
    friend std::ostream& operator<<(std::ostream& out, const BinSettings& b);
};



class AxisSettings : public BinSettings {
protected:
    std::string label;
public:
    AxisSettings(const std::string& label_, const BinSettings& bin_settings) :
        BinSettings(bin_settings), label(label_) {}

    const std::string& Label() const { return label; }
          std::string& Label()       { return label; }
};



class VarBinSettings : public interval<double>  {
protected:
    unsigned int bins;
    std::vector<double> bin_edges;

public:
    struct Exception : std::runtime_error {
        using std::runtime_error::runtime_error;
    };

    /**
     * @brief BinSettings equivalent used for variable bin widths
     * @param edges vector or list containing the bin edges
     */
    VarBinSettings(std::vector<double> edges) :
        interval<double>({edges.front(), edges.back()}),
        bins(edges.size()-1),
        bin_edges(edges)
    {
        if (edges.size() < 2)
            throw Exception("Provided bin edges should contain at least two values");

        if (!IsSane())
            throw Exception("BinSettings: Max < Min");
    }

    unsigned int Bins() const noexcept { return bins; }
    unsigned int& Bins()   noexcept { return bins; }

    std::vector<double> Edges() const noexcept { return bin_edges; }
    std::vector<double>& Edges() noexcept { return bin_edges; }

    /**
     * @brief Get the number of the bin a value would end up in
     * @param v The Value
     * @return a bin number between 0 (lowest bin) to BinSettings::Bins-1, -1 if v is outside the interval
     * @note Bin numbering is different from ROOTs TH{1..3} numbering, where the lowest bin has the index 1
     */
    int getBin(const double v) const noexcept;

    friend std::ostream& operator<<(std::ostream& out, const BinSettings& b);
};


class VarAxisSettings : public VarBinSettings {
protected:
    std::string label;
public:
    /**
     * @brief AxisSettings equivalent to manage axes with a variable bin width
     * @param label_ Axis label
     * @param bin_settings VarBinSettings for the variable bin widths for this axis
     */

    VarAxisSettings(const std::string& label_, const VarBinSettings& bin_settings) :
        VarBinSettings(bin_settings), label(label_) {}
    /**
     * @brief AxisSettings equivalent to manage axes with a variable bin width
     * @param label_ Axis label
     * @param edges list containing the bin edges for the variable bin widths for this axis
     */
    VarAxisSettings(const std::string& label_, const std::initializer_list<double> edges) :
        VarBinSettings(edges), label(label_) {}

    const std::string& Label() const { return label; }
          std::string& Label()       { return label; }
};

}
