/**
  * @brief Test for ProtonPermutation.
  *   Use Candidate cluster sizes to track candidates and particles
  */

#include "catch.hpp"
#include "catch_config.h"

#include <vector>
#include <iostream>
#include <cassert>
#include "analysis/utils/matcher.h"


using namespace std;
using namespace ant;
using namespace ant::analysis;

ostream& operator<<(ostream& stream, const utils::matchpair& m) {

    stream << m.a << "\t->\t" << m.b << ":\t" << m.dist << "\t";

    if(m.matched)
        stream << "matched";
    else
        stream << "not matched";

    return stream;
}

void test_matcher2(const vector<int>& va, const vector<int>& vb, const vector<pair<int,bool>>& exptected) {

    assert(va.size() == exptected.size());

    cout << endl;

    const auto res = utils::match2(va, vb, [] (const int a, const int b) { return abs(a - b);});


    for(const auto& m : res) {
        cout << m << endl;
    }

    for(size_t i=0; i<res.size(); ++i) {
        REQUIRE(res.at(i).a == i);
        REQUIRE(res.at(i).matched == exptected.at(i).second);
        REQUIRE(res.at(i).b == exptected.at(i).first);

    }

    REQUIRE(res.size() == va.size());
}


TEST_CASE("Matcher2: va > vb", "[analysis]") {
    const vector<int> va       = { 30, 40, 10, 20 };
    const vector<int> vb       = { 11, 22, 34 };
    const vector<pair<int,bool>> exp = { {2, true}, {2, false}, {0, true}, {1, true} };

    test_matcher2(va, vb, exp);
}

TEST_CASE("Matcher2: va = vb, all matched", "[analysis]") {
    const vector<int> va       = { 30, 40, 10, 20 };
    const vector<int> vb       = { 11, 22, 34, 39};
    const vector<pair<int,bool>> exp = { {2, true}, {3, true}, {0, true}, {1, true} };

    test_matcher2(va, vb, exp);
}

TEST_CASE("Matcher2: va < vb", "[analysis]") {
    const vector<int> va       = { 30, 40, 10 };
    const vector<int> vb       = { 11, 22, 34, 39};
    const vector<pair<int,bool>> exp = { {2, true}, {3, true}, {0, true} };

    test_matcher2(va, vb, exp);
}
