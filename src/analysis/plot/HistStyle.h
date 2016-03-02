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


// returns string to be used as draw option
struct Mod_t : std::function<std::string(TH1*)> {
    // use constructors of base class
    using std::function<std::string(TH1*)>::function;

    // default ctor modifies nothing
    Mod_t() : std::function<std::string(TH1*)>([] (TH1*) { return ""; }) {}

    // helpers to create modifiers
    static Mod_t MakeLine(Color_t color, short linewidth = 1.0);
    static Mod_t MakeDataPoints(Color_t color, short linewidth = 1.0);

};

}}}}