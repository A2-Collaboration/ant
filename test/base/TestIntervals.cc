#include "catch.hpp"
#include "base/piecewise_interval.h"
#include "base/interval.h"
#include <iostream>

using namespace std;
using namespace ant;

TEST_CASE("Piecewiese Interval: Default ctor", "[base]") {
    REQUIRE_NOTHROW( PiecewiseInterval<int> a; );
}

TEST_CASE("Piecewiese Interval: ctor from interval initializer list", "[base]") {
    REQUIRE_NOTHROW( PiecewiseInterval<int> a({ interval<int>(2,3) }); );
}


TEST_CASE("Piecewiese Interval: Print ", "[base]") {
                PiecewiseInterval<int> a({interval<int>(2,3),interval<int>(5,7)});
                cout << a << endl;
}
TEST_CASE("Piecewiese Interval: Contrains ", "[base]") {
                PiecewiseInterval<int> a({interval<int>(1,3),interval<int>(5,7)});
                REQUIRE(a.Contains(2));
                REQUIRE(!a.Contains(0));
                REQUIRE(a.Contains(6));
                REQUIRE(!a.Contains(9));
}


TEST_CASE("Piecewiese Interval: Compact1 ", "[base]") {
                PiecewiseInterval<int> a({interval<int>(2,4),interval<int>(3,7)});
                cout << a << endl;
                a.Compact();
                cout << a << endl;
                REQUIRE(a.size()==1);
}

TEST_CASE("Piecewiese Interval: Compact2 ", "[base]") {
                PiecewiseInterval<int> a({interval<int>(2,4), interval<int>(0,1), interval<int>(3,7), interval<int>(-5,0)});
                a.Compact();
                REQUIRE(a.size()==2);
                REQUIRE(a == PiecewiseInterval<int>({interval<int>(2,7),interval<int>(-5,1)}));
}
