#pragma once
#include <Rtypes.h>

class TH1;
class TH1D;
class TH2;
class TH2D;
class TF1;
class TCanvas;

namespace ant {
class hstack;

class Pi0 {
public:
    static TCanvas* plotSigmaTheta(const bool showRel = false);
    static TCanvas* plotSigmaThetaMC(const std::string& plotFileName, const bool showRel = false);

    static void plotComparison();


};
}
