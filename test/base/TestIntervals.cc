#include "catch.hpp"
#include "base/piecewise_interval.h"
#include "base/interval.h"
#include "base/interval_algo.h"
#include "base/std_ext/math.h"
#include <iostream>
#include <sstream>

using namespace std;
using namespace ant;

TEST_CASE("Interval: Default ctor", "[base]") {
    REQUIRE_NOTHROW( interval<int> a(0,10); );
}

TEST_CASE("Interval: Length", "[base]") {

    interval<double> a(0,10);
    REQUIRE(a.Length()==10);
}

TEST_CASE("Interval: Inersect", "[base]") {
    interval<int> a(0,100);
    interval<int> b(50,150);
    interval<int> c = intersect(a,b);
    REQUIRE(c.Start()==50);
    REQUIRE(c.Stop()==100);
    interval<int> d = intersect(b,a);
    REQUIRE(c.Start()==50);
    REQUIRE(c.Stop()==100);
}

TEST_CASE("Interval: Parse from string", "[base]") {
    interval<double> a(std_ext::NaN, std_ext::NaN);
    stringstream ss_a;
    ss_a << "[7.3:6.9]";
    REQUIRE(ss_a >> a);
    REQUIRE(a.Start() == Approx(7.3));
    REQUIRE(a.Stop() == Approx(6.9));

    stringstream ss_b;
    ss_b << a;
    interval<double> b(std_ext::NaN, std_ext::NaN);
    REQUIRE(ss_b >> b);
    REQUIRE(b == a);

    stringstream ss_c;
    ss_c << "5.6:2.1 [4.3:6.9]";
    interval<double> c1(std_ext::NaN, std_ext::NaN);
    REQUIRE(ss_c >> c1);
    REQUIRE(c1.Start() == Approx(5.6));
    REQUIRE(c1.Stop() == Approx(2.1));
    interval<double> c2(std_ext::NaN, std_ext::NaN);
    REQUIRE(ss_c >> c2);
    REQUIRE(c2.Start() == Approx(4.3));
    REQUIRE(c2.Stop() == Approx(6.9));

    stringstream ss_d;
    ss_d << "456";
    interval<int> d1(0,0);
    REQUIRE(ss_d >> d1);
    REQUIRE(d1.Start() == 456);
    REQUIRE(d1.Stop() == 456);

    stringstream ss_e;
    ss_e << "6.7:23.1";
    interval<double> e(std_ext::NaN, std_ext::NaN);
    REQUIRE(ss_e >> e);
    REQUIRE(e.Start() == Approx(6.7));
    REQUIRE(e.Stop() == Approx(23.1));
}


TEST_CASE("Piecewiese Interval: Default ctor", "[base]") {
    REQUIRE_NOTHROW( PiecewiseInterval<int> a; );
}

TEST_CASE("Piecewiese Interval: ctor from interval initializer list", "[base]") {
    REQUIRE_NOTHROW( PiecewiseInterval<int> a({ interval<int>(2,3) }); );
}

TEST_CASE("Piecewiese Interval: Print ", "[base]") {
    PiecewiseInterval<int> a({interval<int>(2,3),interval<int>(5,7)});
    stringstream ss_out;
    ss_out << a;
    REQUIRE(ss_out.str() == "[[2:3][5:7]]");
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
    a.Compact();
    REQUIRE(a.size()==1);
    REQUIRE(a.front() == interval<int>(2,7));
}

TEST_CASE("Piecewiese Interval: Compact2 ", "[base]") {
    PiecewiseInterval<int> a({interval<int>(2,4), interval<int>(0,1), interval<int>(3,7), interval<int>(-5,0)});
    a.Compact();
    REQUIRE(a.size()==2);
    REQUIRE(a == PiecewiseInterval<int>({interval<int>(2,7),interval<int>(-5,1)}));
}

TEST_CASE("Piecewiese Interval: Parse from string", "[base]") {
    stringstream ss_a;
    ss_a << "[[2:3][5:7]]";
    PiecewiseInterval<int> a;
    REQUIRE(ss_a >> a);
    REQUIRE(a.size() == 2);
    REQUIRE(a.front() == interval<int>(2,3));
    REQUIRE(a.back() == interval<int>(5,7));

    stringstream ss_b;
    ss_b << "2;6-7";
    PiecewiseInterval<int> b;
    REQUIRE(ss_b >> b);
    REQUIRE(b.size() == 2);
    REQUIRE(b.front() == interval<int>(2,2));
    REQUIRE(b.back() == interval<int>(6,7));

    PiecewiseInterval<int> c_in({{4,5},{3,1},{7,8}});
    stringstream ss_c;
    ss_c << c_in;
    PiecewiseInterval<int> c_out;
    REQUIRE(ss_c >> c_out);
    REQUIRE(c_in == c_out);

}

TEST_CASE("Interval algos", "[base]") {
    const interval<double> a = {0, 10};

    REQUIRE(step(a,3,1) == Approx(5.0));

    std::vector<double> x;
    for(auto s =stepper(a,3); !s.Done(); s.Next()) {
        x.push_back(s.value);
    }

    REQUIRE(x.size() == 3);
    REQUIRE(x.at(0) == Approx(0));
    REQUIRE(x.at(1) == Approx(5));
    REQUIRE(x.at(2) == Approx(10));

    x.clear();

    for(const auto& v : Range(a,3)) {
        x.push_back(v);
    }
    REQUIRE(x.size() == 3);
    REQUIRE(x.at(0) == Approx(0));
    REQUIRE(x.at(1) == Approx(5));
    REQUIRE(x.at(2) == Approx(10));

}
