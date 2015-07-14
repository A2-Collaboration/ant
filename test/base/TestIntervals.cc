#include "catch.hpp"
#include "base/piecewise_interval.h"
#include "base/interval.h"
#include <iostream>

using namespace std;
using namespace ant;

TEST_CASE("Interval: Default ctor", "[base]") {
    REQUIRE_NOTHROW( interval<int> a; );
    REQUIRE_NOTHROW( interval<int> a(0,10); );
}

TEST_CASE("Interval: Length", "[base]") {

    interval<double> a(0,10);
    REQUIRE(a.Length()==10);
}

TEST_CASE("Interval:Inersect", "[base]") {
    interval<int> a(0,100);
    interval<int> b(50,150);
    interval<int> c = intersect(a,b);
    REQUIRE(c.Start()==50);
    REQUIRE(c.Stop()==100);
    interval<int> d = intersect(b,a);
    REQUIRE(c.Start()==50);
    REQUIRE(c.Stop()==100);
}


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
