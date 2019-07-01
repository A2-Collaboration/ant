#include "catch.hpp"

#include "analysis/plot/HistogramFactory.h"
#include "base/WrapTFile.h"
#include "base/tmpfile_t.h"

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
void dotest_numdir();


TEST_CASE("HistogramFactory: Make", "[analysis]") {
    dotest_make();
}


TEST_CASE("HistogramFactory: Name clash", "[analysis]") {
    dotest_nameclash();
}

TEST_CASE("HistogramFactory: Numbered directories", "[analysis]") {
    dotest_numdir();
}


void dotest_make() {
    gDirectory->Clear();

    HistogramFactory h("Test");

    REQUIRE(h.makeTH1D("h",{"",{5,{1,10}}}));
    REQUIRE(h.makeTH1D("h1","","",BinSettings(1)));
    REQUIRE(h.makeTH2D("h2","","",BinSettings(1),BinSettings(2)));
    REQUIRE(h.makeTH3D("h3","","","",BinSettings(1),BinSettings(2),BinSettings(3)));
    REQUIRE(h.makeTH1D("h1_var1", "", "", VarBinSettings({0,1,2})));
    REQUIRE(h.makeTH1D("h1_var2", {"", {0,1,2,3}}));  // c'tor VarAxisSettings
    REQUIRE(h.makeTH1D("h1_var3", {0,1,2,3}, "", ""));  // c'tor vector of edges
    REQUIRE(h.makeTH2D("h2_varx1", "", "", VarBinSettings({0,1,2}), BinSettings(2)));
    REQUIRE(h.makeTH2D("h2_varx2", {"",VarBinSettings({0,1,2})}, {"",BinSettings(2)}));
    REQUIRE(h.makeTH2D("h2_vary1", "", "", BinSettings(2), VarBinSettings({0,1,2})));
    REQUIRE(h.makeTH2D("h2_vary2", {"",2}, {"",{0,1,2,3}}));  // c'tor with AxisSettings and VarAxisSettings
    REQUIRE(h.makeTH2D("h2_varxy1", "", "", VarBinSettings({1,2,3}), VarBinSettings({0,1,2,3})));
    REQUIRE(h.makeTH2D("h2_varxy2", {"x",{0,1,2,3}}, {"y",{1,2,3,4}}));  // c'tor with both axes VarAxisSettings
    REQUIRE(h.makeTTree("tree"));
    REQUIRE(h.makeGraph(""));
    REQUIRE(h.makeGraphErrors(""));
    TH1* h_manual = new TH1D("manual", "", 2, 0, 1);
    REQUIRE_NOTHROW(h.addHistogram(h_manual));
}

void dotest_nameclash() {
    gDirectory->Clear();

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



void dotest_numdir() {
    gDirectory->Clear();

    HistogramFactory h1("Test");
    HistogramFactory h2("Test");
    HistogramFactory h3("Test");

    REQUIRE(dynamic_cast<TDirectory*>(gDirectory->FindObject("Test")));
    REQUIRE(dynamic_cast<TDirectory*>(gDirectory->FindObject("Test_1")));
    REQUIRE(dynamic_cast<TDirectory*>(gDirectory->FindObject("Test_2")));
    REQUIRE(dynamic_cast<TDirectory*>(gDirectory->FindObject("Test_3"))==nullptr);

    // create sub-dir
    HistogramFactory h1_1("SubTest", h1);
    HistogramFactory h1_2("SubTest", h1);
    {
        HistogramFactory::DirStackPush d(h1);
        REQUIRE(dynamic_cast<TDirectory*>(gDirectory->FindObject("SubTest")));
        REQUIRE(dynamic_cast<TDirectory*>(gDirectory->FindObject("SubTest_1")));
        REQUIRE(dynamic_cast<TDirectory*>(gDirectory->FindObject("Test"))==nullptr);
    }

    // back in old dir
    REQUIRE(dynamic_cast<TDirectory*>(gDirectory->FindObject("Test_2")));
}
