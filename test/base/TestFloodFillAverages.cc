#include "catch.hpp"

#include "base/FloodFillAverages.h"
#include "base/std_ext/math.h"

using namespace std;
using namespace ant;
using namespace ant::std_ext;


vector<int> getNeighbours(int i) {
    const vector<vector<int>> neighbours{
        {1,3},   {0,4,2},   {1,5},
        {0,4,6}, {1,3,5,7}, {2,4,8},
        {3,7},   {4,6,8},   {5,7}
    };
    return neighbours[i];
}

void doFloodFill(vector<double>& numbers) {
    floodFillAverages(numbers.size(),
      [&numbers] (int i) { return numbers[i]; },
      [&numbers] (int i, double v) { numbers[i] = v; },
      getNeighbours,
      [&numbers] (int i) { return isfinite(numbers[i]); }
    );
}

void dotest_simple();
void dotest_edge1();
void dotest_edge2();
void dotest_cyclic();


//TEST_CASE("FloodFillAverages: Simple", "[base]") {
//    dotest_simple();
//}

//TEST_CASE("FloodFillAverages: Edge1", "[base]") {
//    dotest_edge1();
//}

//TEST_CASE("FloodFillAverages: Edge2", "[base]") {
//    dotest_edge2();
//}

TEST_CASE("FloodFillAverages: Cyclic", "[base]") {
    dotest_cyclic();
}

void dotest_simple() {
    vector<double> numbers{
        0.5, NaN, 0.5,  // 0,1,2
        0.5, 0.5, 0.5,  // 3,4,5
        NaN, NaN, 0.5   // 6,7,8
    };

    doFloodFill(numbers);

    for(auto& n : numbers)
        CHECK(n == Approx(0.5));
}


void dotest_edge1() {
    vector<double> numbers{
        NaN, NaN, NaN,  // 0,1,2
        NaN, NaN, NaN,  // 3,4,5
        NaN, NaN, 0.5   // 6,7,8
    };

    doFloodFill(numbers);

    for(auto& n : numbers)
        CHECK(n == Approx(0.5));

}

void dotest_edge2() {
    vector<double> numbers{
        0.1, NaN, NaN,  // 0,1,2
        NaN, NaN, NaN,  // 3,4,5
        NaN, NaN, 0.5   // 6,7,8
    };

    doFloodFill(numbers);

    CHECK(numbers[1] == Approx(0.1));
    CHECK(numbers[3] == Approx(0.1));
    CHECK(numbers[2] == Approx(0.3));
    CHECK(numbers[4] == Approx(0.3));
    CHECK(numbers[6] == Approx(0.3));
    CHECK(numbers[7] == Approx(0.5));
    CHECK(numbers[5] == Approx(0.5));
}

void dotest_cyclic() {
    vector<double> ring{
        0.1, NaN, NaN, NaN, NaN, 0.2, NaN
    };

    floodFillAverages(ring.size(),
      [&ring] (int i) { return ring[i]; },
      [&ring] (int i, double v) { ring[i] = v; },
      [&ring] (int i) -> vector<int> { if(i==0) return {}; return {i-1};},
      [&ring] (int i) { return isfinite(ring[i]); }
    );

    for(unsigned i=0;i<ring.size();i++)
        CHECK(ring[i] == Approx(i<5 ? 0.1 : 0.2));

}
