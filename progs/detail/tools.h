#pragma once

#include "base/piecewise_interval.h"
#include "base/std_ext/string.h"

#include <stdexcept>

namespace ant {
namespace progs {
namespace tools {

PiecewiseInterval<unsigned> parse_cmdline_ranges(const std::vector<std::string>& strings) {
    PiecewiseInterval<unsigned> ranges;
    for(auto& str : strings) {
        std::stringstream ss(str);
        interval<unsigned> range(0,0);
        if(ss >> range && range.IsSane())
            ranges.push_back(range);
        else
            throw std::runtime_error(std_ext::formatter() << "Cannot parse cmdline input range '" << str << "'");
    }
    return ranges;
};

}}} // namespace ant::progs::tools
