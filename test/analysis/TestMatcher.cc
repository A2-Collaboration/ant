/**
  * @brief Test for ProtonPermutation.
  *   Use Candidate cluster sizes to track candidates and particles
  */

#include "catch.hpp"
#include "catch_config.h"

#include <vector>
#include <iostream>
#include "analysis/utils/matcher.h"


using namespace std;
using namespace ant;
using namespace ant::analysis;

void test_matcher2();

TEST_CASE("Analysis: Matcher2", "[analysis]") {
    test_matcher2();
}

ostream& operator<<(ostream& stream, const utils::matchpair& m) {

    stream << m.a << "\t->\t" << m.b << ":\t" << m.dist << "\t";

    if(m.matched)
        stream << "matched";
    else
        stream << "not matched";

    return stream;
}

void test_matcher2() {
    const vector<int> va = { 30, 40, 10, 20 };
    const vector<int> vb = { 11, 22, 34 };

    // 40 will be unmatched, but closest to 34.

    const auto res = utils::match2(va, vb, [] (const int a, const int b) { return abs(a - b);});

//    for(const auto& m : res) {
//        cout << m << endl;
//    }

    REQUIRE(res.at(0).matched == true);
    REQUIRE(res.at(0).b == 2);
    REQUIRE(res.at(1).matched == false);
    REQUIRE(res.at(1).b == 2);
    REQUIRE(res.at(2).matched == true);
    REQUIRE(res.at(2).b == 0);
    REQUIRE(res.at(3).matched == true);
    REQUIRE(res.at(3).b == 1);


}
