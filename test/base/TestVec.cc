#include "catch.hpp"

#include "base/vec3.h"
#include "base/LorentzVec.h"

#include <type_traits>
#include <iostream>


using namespace std;
using namespace ant;

void do_test_3();
void do_test_lv();

TEST_CASE("vec3", "[base]") {
    do_test_3();
}

TEST_CASE("LorentzVec", "[base]") {
    do_test_lv();
}

void do_test_3() {
    REQUIRE(std::is_constructible<vec3>::value);
    REQUIRE(std::is_default_constructible<vec3>::value);
    REQUIRE(std::is_nothrow_destructible<vec3>::value);
    REQUIRE(std::is_copy_constructible<vec3>::value);
    REQUIRE(std::is_move_constructible<vec3>::value);
    REQUIRE(std::is_nothrow_copy_assignable<vec3>::value);
    REQUIRE(std::is_nothrow_move_assignable<vec3>::value);

    const TVector3 a(1,2,3);
    const vec3     b(1,2,3);

    REQUIRE(a == b);

    const auto c = b;
    REQUIRE(b==c);
    REQUIRE(a==c);

    REQUIRE(a.Mag()   == Approx(b.R()));
    REQUIRE(a.Theta() == Approx(b.Theta()));
    REQUIRE(a.Phi()   == Approx(b.Phi()));

    const TVector3 f = b;
    REQUIRE(a == f);

    REQUIRE( b*3 == vec3(3,6,9));
    REQUIRE( 3*b == vec3(3,6,9));

    auto x = b;
    x*=3;
    REQUIRE(x == vec3(3,6,9));
    x/=3;
    REQUIRE(x == b);

    auto y = x+b;
    REQUIRE(y == vec3(2,4,6));

    REQUIRE(b.Dot(y)   == Approx(a.Dot(y)));
    REQUIRE(b.Angle(y) == Approx(a.Angle(y)));

    auto m = vec3::RThetaPhi(1.0, 2.0, 3.0);
    REQUIRE(m.R()     == Approx(1.0));
    REQUIRE(m.Theta() == Approx(2.0));
    REQUIRE(m.Phi()   == Approx(3.0));
    m *= 5.0;
    REQUIRE(m.R()     == Approx(5.0));
    REQUIRE(m.Theta() == Approx(2.0));
    REQUIRE(m.Phi()   == Approx(3.0));

}

void do_test_lv() {
    REQUIRE(std::is_constructible<LorentzVec>::value);
    REQUIRE(std::is_default_constructible<LorentzVec>::value);
    REQUIRE(std::is_nothrow_destructible<LorentzVec>::value);
    REQUIRE(std::is_copy_constructible<LorentzVec>::value);
    REQUIRE(std::is_move_constructible<LorentzVec>::value);
    REQUIRE(std::is_nothrow_copy_assignable<LorentzVec>::value);
    REQUIRE(std::is_nothrow_move_assignable<LorentzVec>::value);

    const TLorentzVector a(1,2,3,4);
    const LorentzVec     b(1,2,3,4);
    const TLorentzVector a2(3,2,4,4);
    const LorentzVec     b2(3,2,4,4);

    REQUIRE(a == b);

    const auto c = b;
    REQUIRE(b==c);
    REQUIRE(a==c);

    REQUIRE(a.M() == Approx(b.M()));
    REQUIRE(a.Theta() == Approx(b.Theta()));
    REQUIRE(a.Phi() == Approx(b.Phi()));
    REQUIRE(a.Beta() == Approx(b.Beta()));
    REQUIRE(a.Gamma() == Approx(b.Gamma()));
    REQUIRE(a.BoostVector() == b.BoostVector());
    REQUIRE(b.Dot(b2) == Approx(a.Dot(a2)));
    REQUIRE(a.P() == Approx(b.P()));

    const auto x = a+a2;
    const auto y = b+b2;

    REQUIRE(a.M() == Approx(b.M()));

    REQUIRE(a+a2 == b+b2);
    REQUIRE(a-a2 == b-b2);
    REQUIRE(a*3.0 == b*3.0);
    REQUIRE(a*(1./3.0) == b/3.0);

    REQUIRE(a.Vect() == b.x);
}
