#include "BinSettings.h"
#include <cmath>
#include <numeric>

using namespace std;
using namespace ant;

BinSettings BinSettings::RoundToBinSize(const BinSettings& bins, const double binSize) {

    unsigned nOptBins = unsigned(bins.Length() / binSize);

    const auto factor = unsigned(std::round(double(nOptBins) / double(bins.Bins())));

    if(nOptBins >= binSize) {
         nOptBins /= factor;
    } else {
        nOptBins *= factor;
    }

    const auto length = binSize * nOptBins * factor;

    return BinSettings(nOptBins, interval<double>::CenterWidth(bins.Center(), length));
}

BinSettings BinSettings::Make(const vector<double>& x_values)
{
    if(x_values.size()<2)
        throw runtime_error("Given x_values should contain at least 2 elements");
    vector<double> diffs(x_values.size());
    std::adjacent_difference(x_values.begin(), x_values.end(), diffs.begin());
    if(std::find_if(diffs.begin(), diffs.end(), [] (double v) { return v<=0; }) != diffs.end())
        throw runtime_error("Given x_values not strictly monotonic");
    double mean_distance = std::accumulate(std::next(diffs.begin()), diffs.end(), 0.0) / diffs.size();
    return { static_cast<unsigned>(x_values.size()),
                x_values.front() - mean_distance/2.0,
                x_values.back()  + mean_distance/2.0 };
}

int BinSettings::getBin(const double v) const noexcept {
    if(this->Contains(v)) {
        const auto x = v - this->Start();
        const auto t = x / this->Length(); // [0..1]
        const auto bin = bins * t; // [0..bins]
        return int(std::floor(bin));
    } else
        return -1;
}


int VarBinSettings::getBin(const double v) const noexcept {
    if (!this->Contains(v))
        return -1;

    size_t pos = 0;
    while (pos < bins) {
        if (bin_edges.at(pos) <= v && bin_edges.at(pos+1) < v)  // current and next edge below value
            pos++;  // increase position to check
        else
            break;  // above if condition evaluated to false: found position where lowEdge <= v < highEdge
    }
    return int(pos);
}



namespace ant {

istream& operator>>(istream& in, BinSettings& t)
{
    // skip leading whitespace
    in >> std::ws;
    // skip leading (
    if(in.peek() == '(') {
        in.ignore();
    }
    // read Start
    if(in >> t.Bins()) {
        // maybe read stop?
        if(in.peek() == ',') {
            in.ignore();
            in >> static_cast<interval<double>&>(t);
        }
        // optionally have closing )
        // do not mark EOF as error
        auto nextchar = in.peek();
        if(nextchar == ')')
            in.ignore();
        else if(nextchar == EOF)
            in.clear();
    }
    return in;
}

ostream& operator<<(ostream& out, const BinSettings& b)
{
    out << "(" << b.Bins() << "," << static_cast<const ant::interval<double>&>(b) << ")";
    return out;
}

ostream& operator<<(ostream& out, const VarBinSettings& b)
{
    out << "variable bins (" << b.Bins() << "," << static_cast<const ant::interval<double>&>(b) << ")";
    return out;
}

} // namespace ant
