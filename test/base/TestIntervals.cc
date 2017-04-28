#include "catch.hpp"
#include "base/piecewise_interval.h"
#include "base/interval.h"
#include "base/interval_algo.h"
#include "base/std_ext/math.h"
#include "base/std_ext/memory.h"
#include <iostream>
#include <sstream>

using namespace std;
using namespace ant;

TEST_CASE("Interval: Default ctor", "[base]") {
    REQUIRE_NOTHROW( std_ext::make_unique<interval<int>>(0,10) );
}

TEST_CASE("Interval: Length", "[base]") {

    interval<double> a(0,10);
    REQUIRE(a.Length()==10);
}

TEST_CASE("Interval: Intersect", "[base]") {
    interval<int> a(0,100);
    interval<int> b(50,150);
    interval<int> c = intersect(a,b);
    REQUIRE(c.Start()==50);
    REQUIRE(c.Stop()==100);
    interval<int> d = intersect(b,a);
    REQUIRE(d.Start()==50);
    REQUIRE(d.Stop()==100);
}

TEST_CASE("Interval: Disjoint", "[base]") {
    interval<int> a(0,10);
    interval<int> b(5,10);
    interval<int> c(-5,20);
    interval<int> d(-5, -3);

    REQUIRE(!a.Disjoint(a));
    REQUIRE(!a.Disjoint(b));
    REQUIRE(!a.Disjoint(c));
    REQUIRE( a.Disjoint(d));

    REQUIRE(!b.Disjoint(a));
    REQUIRE(!b.Disjoint(b));
    REQUIRE(!b.Disjoint(c));
    REQUIRE( b.Disjoint(d));

    REQUIRE(!c.Disjoint(a));
    REQUIRE(!c.Disjoint(b));
    REQUIRE(!c.Disjoint(c));
    REQUIRE(!c.Disjoint(d));

    REQUIRE( d.Disjoint(a));
    REQUIRE( d.Disjoint(b));
    REQUIRE(!d.Disjoint(c));
    REQUIRE(!d.Disjoint(d));
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

TEST_CASE("Interval: Parse from string fails", "[base]") {
    interval<double> a(std_ext::NaN, std_ext::NaN);

    {
        stringstream ss;
        ss << "[7.3,6.9]";
        REQUIRE_FALSE(ss >> a);
    }

    {
        stringstream ss;
        ss << "[7.3:6,9]";
        REQUIRE_FALSE(ss >> a);
    }

    {
        stringstream ss;
        ss << "7.3,6.9";
        REQUIRE_FALSE(ss >> a);
    }

    {
        stringstream ss;
        ss << "[7.3:6.9";
        REQUIRE_FALSE(ss >> a);
    }

    {
        stringstream ss;
        ss << "7.3:6.9]";
        REQUIRE_FALSE(ss >> a);
    }
}

TEST_CASE("Interval Clip", "[base]") {
    constexpr interval<double> i = {1,2};
    REQUIRE(i.Clip(1) == 1);
    REQUIRE(i.Clip(1.5) == 1.5);
    REQUIRE(i.Clip(-1) == 1);
    REQUIRE(i.Clip(2) == 2);
    REQUIRE(i.Clip(3) == 2);
    REQUIRE(i.Clip(std_ext::inf) == 2);
    REQUIRE(i.Clip(-std_ext::inf) == 1);
    REQUIRE(std::isnan(i.Clip(-std_ext::NaN)));
    REQUIRE(std::isnan(i.Clip(+std_ext::NaN)));
}

TEST_CASE("Interval AsCutString", "[base]") {
    REQUIRE(interval<double>(1,2).AsRangeString() == "1<=x<=2");
    REQUIRE(interval<double>(1,2).AsRangeString("y") == "1<=y<=2");
    REQUIRE(interval<double>(std_ext::NaN,2).AsRangeString("y") == "y<=2");
    REQUIRE(interval<double>(-std_ext::inf,2).AsRangeString("y") == "y<=2");
    REQUIRE(interval<double>(1,std_ext::NaN).AsRangeString("y") == "1<=y");
    REQUIRE(interval<double>(1,std_ext::inf).AsRangeString("y") == "1<=y");
    REQUIRE(interval<double>(-std_ext::inf,std_ext::inf).AsRangeString("y") == "y");
    REQUIRE(interval<double>(std_ext::inf,-std_ext::inf).AsRangeString("y") == "y");
}

TEST_CASE("Interval Round", "[base]") {
    REQUIRE(IntervalD(1.2,1.0).Round() == IntervalD(1.0, 1.0));
    REQUIRE(IntervalD(1.2,2.2).Round() == IntervalD(1.0, 2.0));
    REQUIRE(IntervalD(-std_ext::inf,1.5).Round() == IntervalD(-std_ext::inf, 2.0));
    REQUIRE(std::isnan(IntervalD(-std_ext::inf,std_ext::NaN).Round().Stop()));
}

TEST_CASE("Piecewiese Interval: Default ctor", "[base]") {
    REQUIRE_NOTHROW( std_ext::make_unique<PiecewiseInterval<int>>() );
}

TEST_CASE("Piecewiese Interval: ctor from interval initializer list", "[base]") {
    const std::initializer_list<interval<int>> i{interval<int>(2,3)};
    REQUIRE_NOTHROW( std_ext::make_unique<PiecewiseInterval<int>>(i) );
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
