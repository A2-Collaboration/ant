#include "catch.hpp"

#include "analysis/plot/HistogramFactory.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TGraph.h"
#include "TTree.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;

void dotest_make();
void dotest_nameclash();

TEST_CASE("HistogramFactory: Make", "[analysis]") {
    dotest_make();
}


TEST_CASE("HistogramFactory: Name clash", "[analysis]") {
    dotest_nameclash();
}

void dotest_make() {
    HistogramFactory h("Test");

    REQUIRE(h.makeTH1D("h1","","",BinSettings(1)));
    REQUIRE(h.makeTH2D("h2","","",BinSettings(1),BinSettings(2)));
    REQUIRE(h.makeTH3D("h3","","","",BinSettings(1),BinSettings(2),BinSettings(3)));
    REQUIRE(h.makeTTree("tree"));
    REQUIRE(h.makeGraph(""));
    REQUIRE(h.makeGraphErrors(""));
}

void dotest_nameclash() {
    HistogramFactory h("Test");
    HistogramFactory h2("Test2", h);
    REQUIRE_NOTHROW(h.makeTTree("t1"));
    REQUIRE_NOTHROW(h2.makeTTree("t1"));
    REQUIRE_THROWS_AS(h.makeTTree("t1"), HistogramFactory::Exception);
    REQUIRE_THROWS_AS(h.makeTTree(""), HistogramFactory::Exception);
    REQUIRE_NOTHROW(h.makeTH1D("h1","","",BinSettings(1)));
    REQUIRE_NOTHROW(h.makeTH1D("h1","","",BinSettings(1)));
    REQUIRE_NOTHROW(h.makeGraph(""));
    REQUIRE_NOTHROW(h.makeGraph(""));
}
