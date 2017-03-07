#include "catch.hpp"

#include "base/SavitzkyGolay.h"

using namespace std;
using namespace ant;


TEST_CASE("SavitzkyGolay: Simple average", "[base/std_ext]") {
    SavitzkyGolay sg(2,2,0);
    auto smoothed = sg.Smooth({1,1,1,1,1});
    // zero-padding, args
    auto expected = vector<double>{0.6,0.8,1,0.8,0.6};
    for(auto i=0u;i<expected.size();i++)
        REQUIRE(smoothed[i] == Approx(expected[i]));
}
