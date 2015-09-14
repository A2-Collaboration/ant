#include "plot/SmartHist.h"
#include "TH1D.h"

using namespace ant::analysis;

template<>
SmartHist1<std::string> SmartHist1<std::string>::makeHist(
    const std::string& title,
    const std::string& xlabel,
    const std::string& ylabel,
    const BinSettings& bins,
    const std::string& name,
    HistogramFactory& factory)
{
    auto converter = [] (const std::string& data) -> const char* {
        return data.c_str();
    };
    return makeHist(converter, title, xlabel, ylabel, bins, name, factory);
}
