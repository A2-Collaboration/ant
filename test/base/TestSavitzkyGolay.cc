#include "catch.hpp"

#include "base/SavitzkyGolay.h"

using namespace std;
using namespace ant;


TEST_CASE("SavitzkyGolay: Simple average", "[base/std_ext]") {
    SavitzkyGolay sg(2,2,0);
    auto smoothed = sg.Smooth({1,1,1,1,1});
    auto expected = vector<double>{1,1,1,1,1};
    REQUIRE(smoothed == expected);
}
