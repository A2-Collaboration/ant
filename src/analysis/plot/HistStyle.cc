#include "HistStyle.h"

#include "base/Logger.h"

#include "TH1.h"


#include <cassert>

using namespace std;
using namespace ant;
using namespace ant::analysis::plot::histstyle;


Color_t color_t::Get(unsigned i) {
    static std::vector<Color_t> colors;
    if(colors.empty()) {
        colors = {
            // have a look at the TColorWheel for this list
            kOrange-9,
            kYellow-9,
            kSpring-8,
            kGreen-9,
            kTeal-9,
            kCyan-9,
            kAzure-9,
            kBlue-9,
            kViolet-9,
            kMagenta-9,
            kPink-9,
            kRed-9
        };
    }
    return colors[i % colors.size()];
}

Mod_t Mod_t::MakeLine(Color_t color, short linewidth) {
    return [color, linewidth] (TH1* h) {
        h->SetLineColor(color);
        h->SetLineWidth(linewidth);
        h->SetMarkerSize(1);
        h->SetMarkerColor(color);
        return "";
    };
}

Mod_t Mod_t::MakeDataPoints(Color_t color, short linewidth) {
    return [color, linewidth] (TH1* h) {
        h->SetLineColor(color);
        h->SetLineWidth(linewidth);
        h->SetMarkerColor(color);
        h->SetMarkerStyle(kDot);
        return "E"; // draw error bars
    };

}
