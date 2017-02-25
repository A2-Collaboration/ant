#pragma once

#include "Rtypes.h" // for Color_t

#include <string>
#include <functional>



class TH1;

namespace ant {
namespace analysis {
namespace plot {
namespace histstyle {

struct color_t {
    /**
     * @brief Get returns well-chosen color for dark/greyish backgrounds
     * @param i
     * @return
     */
    static Color_t GetLight(unsigned i);
    /**
     * @brief GetDark returns well-chosen color for bright/whitish backgrounds
     * @param i
     * @return
     */
    static Color_t GetDark(unsigned i);
};


struct ModOption_t {
    std::string DrawOption;
    int Z;
    Color_t BkgColor;
    ModOption_t(const std::string& drawoption = "", int z = 0, Color_t bkgColor = -1) :
        DrawOption(drawoption), Z(z), BkgColor(bkgColor)
    {}
    template<typename Archive>
    void serialize(Archive archive) {
        archive(DrawOption, Z, BkgColor);
    }
};

struct Mod_t : std::function<ModOption_t(TH1*)> {
    // use constructors of base class
    using std::function<ModOption_t(TH1*)>::function;

    // default ctor modifies nothing
    Mod_t() : std::function<ModOption_t(TH1*)>([] (TH1*) { return ModOption_t{}; }) {}

    // helpers to create modifiers
    static Mod_t MakeLine(Color_t color, short linewidth = 1, Color_t bkgColor = -1);
    static Mod_t MakeDataPoints(Color_t color, short linewidth = 1);
    static Mod_t MakeFill(Color_t color, int zpos = -1);

};

}}}}