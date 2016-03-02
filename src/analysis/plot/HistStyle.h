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
    static Color_t Get(unsigned i);
};


struct ModOption_t {
    std::string DrawOption;
    int Z;
    ModOption_t(const std::string& drawoption = "", int z = 0) :
        DrawOption(drawoption), Z(z)
    {}
    template<typename Archive>
    void serialize(Archive archive) {
        archive(DrawOption, Z);
    }
};

struct Mod_t : std::function<ModOption_t(TH1*)> {
    // use constructors of base class
    using std::function<ModOption_t(TH1*)>::function;

    // default ctor modifies nothing
    Mod_t() : std::function<ModOption_t(TH1*)>([] (TH1*) { return ModOption_t{}; }) {}

    // helpers to create modifiers
    static Mod_t MakeLine(Color_t color, short linewidth = 1);
    static Mod_t MakeDataPoints(Color_t color, short linewidth = 1);
    static Mod_t MakeFill(Color_t color, int zpos = -1);

};

}}}}