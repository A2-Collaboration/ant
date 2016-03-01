#pragma once

#include "base/std_ext/math.h"

#include "TH1.h"
#include "TColor.h"

#include <vector>
#include <functional>

namespace ant {
namespace analysis {
namespace plot {
namespace histstyle {

struct color_t {
    double R = std_ext::NaN;
    double G = std_ext::NaN;
    double B = std_ext::NaN;

    color_t(double r, double g, double b);
    color_t(const Color_t color);

    Color_t ToColor_t() const;

    static color_t Get(unsigned i);

protected:
    // count instances
    static unsigned n;
    Color_t root_color;
};


// returns string to be used as draw option
struct Mod_t : std::function<std::string(TH1*)> {
    // use constructors of base class
    using std::function<std::string(TH1*)>::function;

    Mod_t() : std::function<std::string(TH1*)>([] (TH1*) { return ""; }) {}

    static Mod_t MakeLine(const color_t& color, short linewidth = 1.0) {
        return [&color, linewidth] (TH1* h) {
            h->SetLineColor(color.ToColor_t());
            h->SetLineWidth(linewidth);
            h->SetMarkerSize(1);
            h->SetMarkerColor(color.ToColor_t());
            return "";
        };
    }

    static Mod_t MakeDataPoints(const color_t& color, short linewidth = 1.0) {
        return [color, linewidth] (TH1* h) {
            h->SetLineColor(color.ToColor_t());
            h->SetLineWidth(linewidth);
            h->SetMarkerColor(color.ToColor_t());
            h->SetMarkerStyle(kDot);
            return "E"; // draw error bars
        };

    }

};

}}}}