#include "catch.hpp"

#include "analysis/plot/root_draw.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TGraph.h"
#include "TTree.h"
#include "TDirectory.h"

using namespace std;
using namespace ant;

void dotest_make();
void dotest_nameclash();

TEST_CASE("TTreeDrawable: Make", "[analysis]") {
    dotest_make();
}

void dotest_make() {

    TTree *t = new TTree("t","default tree");
    int x;
    int y;

    t->Branch("x",&x);
    t->Branch("y",&y);

    for (auto i = 0; i < 10; i++){
        x = i;
        y = i;
        t->Fill();
    }

    // Test creation of basic histograms
    // ============================================================
    delete gDirectory->Get("htemp0");        // Delete just in case
    delete gDirectory->Get("htemp1");        // Delete just in case
    ant::canvas("Test")
            << TTree_drawable(t, "x", "")       // Should be htemp0
            << TTree_drawable(t, "y", "")       // Should be htemp1
            << endc;

    TH1* h0; gDirectory->GetObject("htemp0",h0);
    TH1* h1; gDirectory->GetObject("htemp1",h1);

    REQUIRE(h0);
    REQUIRE(h1);

    // Test setting bins
    // ============================================================
    delete gDirectory->Get("test_bins");
    ant::canvas("Test")
            << TTree_drawable(t, "x>>test_bins(20,0,20)", "")
            << endc;

    TH1* h2; gDirectory->GetObject("test_bins",h2);
    REQUIRE(h2);
    REQUIRE(h2->GetNbinsX() == 20);


    // Test advanced constructor
    // ============================================================
    // TTree_drawable(TTree* tree, const std::string& varx, const std::string& cut,  const std::string &title, const std::string &xlabel, const std::string &ylabel, const BinSettings &binsx);

    delete gDirectory->Get("test_const");
    ant::canvas("Test")
            << TTree_drawable(t, "x", "", "test_const","testx","testy",BinSettings(20))
            << endc;

    TH1* h3; gDirectory->GetObject("test_const",h3);
    REQUIRE(h3);
    REQUIRE(h3->GetNbinsX() == 20);
    REQUIRE(strcmp(h3->GetXaxis()->GetTitle(),"testx") == 0);
    REQUIRE(strcmp(h3->GetYaxis()->GetTitle(),"testy") == 0);

}
