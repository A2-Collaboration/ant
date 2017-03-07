#include "catch.hpp"

#include "base/SavitzkyGolay.h"

using namespace std;
using namespace ant;


TEST_CASE("SavitzkyGolay: Constant values", "[base/std_ext]") {
    SavitzkyGolay sg(2,2,3); // needs at least 2+2+1 points
    auto input = vector<double>{1, 1, 1, 1, 1};
    auto smoothed = sg.Smooth(input);
    auto expected = input;
    for(auto i=0u;i<expected.size();i++) {
        INFO(i);
        REQUIRE(smoothed[i] == Approx(expected[i]));
    }
}
