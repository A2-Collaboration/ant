#include "catch.hpp"
#include "base/bitflag.h"

using namespace std;
using namespace ant;

TEST_CASE("bitflag", "[base]") {

    enum class flags_t {
        F1, F2, F3, F4
    };

    using myflags_t = bitflag<flags_t>;

    myflags_t f1;
    myflags_t f2;

    REQUIRE_FALSE((f1 & flags_t::F1));
    REQUIRE_FALSE((flags_t::F1 & f1));

}
