#include "plot/SmartHist.h"
#include "TH1D.h"

template<>
ant::SmartHist1<std::string> ant::SmartHist1<std::string>::makeHist(
    const std::string& title,
    const std::string& xlabel,
    const std::string& ylabel,
    const BinSettings& bins,
    const std::string& name,
    HistogramFactory& factory)
{
    return move(makeHist([] (const std::string& data) -> const char* { return data.c_str();}, title, xlabel, ylabel, bins, name, factory));
}
