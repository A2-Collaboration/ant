#include "catch.hpp"
#include "catch_config.h"

#include "Unpacker.h"
#include "expconfig/ExpConfig.h"

#include "tree/TEvent.h"

#include "base/tmpfile_t.h"
#include "base/std_ext/misc.h"

#include "TFile.h"
#include "TTree.h"

#include <iostream>
#include <string>

// some REQUIRE statements produce those warnings
// ignore them in this file
#pragma GCC diagnostic ignored "-Wparentheses"

using namespace std;
using namespace ant;

void dotest();

TEST_CASE("Test TreeWriter TEvent's", "[unpacker]") {
    dotest();
}

void dotest() {

    /// \todo the ROOT object handling here is baaaad,
    /// the test segfaults when some requirements are not met

    ant::tmpfile_t tmpfile;

    ant::ExpConfig::Setup::ManualName = "Setup_Test";
    auto unpacker = ant::Unpacker::Get(string(TEST_BLOBS_DIRECTORY)+"/Acqu_oneevent-big.dat.xz");

    // write some stuff to a ROOT tree
    TFile* file = new TFile(tmpfile.filename.c_str(),"RECREATE");

    TTree* tree = new TTree("tree","tree");
    TEvent* ptr = new TEvent();
    REQUIRE(tree->Branch("branch", "ant::TEvent", &ptr) != nullptr);

    while(auto item = unpacker->NextEvent()) {
        ptr = item.get();
        tree->Fill();
    }

    REQUIRE( tree->GetEntries() == 221 );

    file->Write();
    file->Close();
    delete file;

    // open it again

    file = new TFile(tmpfile.filename.c_str());

    REQUIRE(file->IsOpen());

    tree = dynamic_cast<TTree*>(file->Get("tree"));
    REQUIRE((tree != nullptr));
    REQUIRE((tree->GetEntries() == 221));

    ptr = 0;
    REQUIRE(tree->SetBranchAddress("branch", &ptr) == 0);

    for(Long64_t i=0;i<tree->GetEntries();i++) {
        delete ptr;
        ptr = 0;
        tree->GetEntry(i);
        //cout << ptr << endl;
    }

    file->Close();

}
