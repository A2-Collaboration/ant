#include "HistStyle.h"

#include "base/Logger.h"

#include "TROOT.h"

#include <cassert>

using namespace std;
using namespace ant;
using namespace ant::analysis::plot::histstyle;

unsigned color_t::n = 0;

color_t color_t::Get(unsigned i) {
    static std::vector<color_t> colors;
    if(colors.empty()) {
        unsigned n = 0;
        auto rgb = [&n] (double r, double g, double b) {
            return color_t{r/255,g/255,b/255};
        };
        colors = {
            // color values copied from http://colorbrewer2.org/
            rgb(141,211,199),
            rgb(255,255,179),
            rgb(190,186,218),
            rgb(251,128,114),
            rgb(128,177,211),
            rgb(253,180,98),
            rgb(179,222,105),
            rgb(252,205,229),
            rgb(217,217,217),
            rgb(188,128,189),
            rgb(204,235,197),
            rgb(255,237,111),
        };
    }
    return colors[i % colors.size()];
}


color_t::color_t(double r, double g, double b) :
    R(r), G(g), B(b)
{
    root_color = 2000 + n++;
}

color_t::color_t(const Color_t color)
{
    root_color = color;
}

Color_t color_t::ToColor_t() const
{
    if(gROOT->GetColor(root_color) == nullptr) {
        LOG(INFO) << "Adding TColor: idx=" << root_color << " " << R<<" "<<G<<" "<<B;
        new TColor(root_color, R, G, B);
    }
    assert(gROOT->GetColor(root_color) != nullptr);
    return root_color;
}
