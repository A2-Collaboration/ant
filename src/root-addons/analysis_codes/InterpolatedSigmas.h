#pragma once
#include <string>

class TH2;
class TDirectory;

namespace ant {

struct InterpolatedSigamas {

    static void PlotDirectory(TDirectory* dir, const std::string& prefix, const int cols, const int rows, const std::string& title="", TDirectory* dir2=nullptr);
    static void PlotPullSigmas(const std::string& treename, TDirectory* dir);

    static void PlotComparePulls(TDirectory* red, TDirectory* blue);

};

}
