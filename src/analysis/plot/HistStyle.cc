#include "HistStyle.h"

#include "base/Logger.h"

#include "TH1.h"


#include <cassert>

using namespace std;
using namespace ant;
using namespace ant::analysis::plot::histstyle;


Color_t color_t::GetLight(unsigned i) {
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

Color_t color_t::GetDark(unsigned i)
{
    static std::vector<Color_t> colors;
    if(colors.empty()) {
        colors = {
            // have a look at the TColorWheel for this list
            kBlue,
            kRed,
            kGreen+1,
            kBlack,
            kMagenta+1,
            kCyan+1,
            kOrange+7,
        };
    }
    return colors[i % colors.size()];
}

Mod_t Mod_t::MakeLine(Color_t color, short linewidth, Color_t bkgColor) {
    return [color, linewidth, bkgColor] (TH1* h) {
        h->SetLineColor(color);
        h->SetLineWidth(linewidth);
        h->SetMarkerSize(1);
        h->SetMarkerColor(color);
        return ModOption_t{"", 0, bkgColor};
    };
}

Mod_t Mod_t::MakeDataPoints(Color_t color, short linewidth) {
    return [color, linewidth] (TH1* h) {
        h->SetLineColor(color);
        h->SetLineWidth(linewidth);
        h->SetMarkerColor(color);
        h->SetMarkerStyle(kDot);
        if(h->GetSumw2N()==0)
            h->Sumw2();
        return ModOption_t{"E"}; // draw error bars
    };

}

/**
 * @brief Mod_t::MakeFill
 * @param color
 * @param zpos: higher -> in front
 * @return
 */
Mod_t Mod_t::MakeFill(Color_t color, int zpos)
{
    return [color, zpos] (TH1* h) {
        h->SetLineWidth(1);
        h->SetLineColor(color);
        h->SetMarkerColor(color);
        h->SetFillColor(color);
        h->SetFillStyle(1001); // using kFSolid is wrong...
        return ModOption_t{"", zpos};
    };
}
