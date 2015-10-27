#include "catch.hpp"
#include "catch_config.h"

#include "Unpacker.h"
#include "expconfig/ExpConfig.h"

#include "tree/THeaderInfo.h"
#include "tree/TUnpackerMessage.h"
#include "tree/TSlowControl.h"
#include "tree/TDetectorRead.h"

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

template<typename T>
void dotest(Long64_t expectedEntries);

TEST_CASE("Test TreeWriter: ant::THeaderInfo", "[unpacker]") {
    dotest<ant::THeaderInfo>(1);
}

TEST_CASE("Test TreeWriter: ant::TSlowControl", "[unpacker]") {
    dotest<ant::TSlowControl>(8);
}

TEST_CASE("Test TreeWriter: ant::TUnpackerMessage", "[unpacker]") {
    dotest<ant::TUnpackerMessage>(8);
}

TEST_CASE("Test TreeWriter: ant::TDetectorRead", "[unpacker]") {
    dotest<ant::TDetectorRead>(221);
}

template<typename T>
void dotest(Long64_t expectedEntries) {

    const string& classname = ant::std_ext::getTypeAsString<T>();

    ant::tmpfile_t tmpfile;

    ant::ExpConfig::Setup::ManualName = "Setup_Test";
    auto unpacker = ant::Unpacker::Get(string(TEST_BLOBS_DIRECTORY)+"/Acqu_oneevent-big.dat.xz");

    // write some stuff to a ROOT tree
    TFile* file = new TFile(tmpfile.filename.c_str(),"RECREATE");

    TTree* tree = new TTree("tree","tree");
    T* ptr = new T();
    REQUIRE(tree->Branch("branch", classname.c_str(), &ptr) != nullptr);

    while(auto item = unpacker->NextItem()) {
        ptr = dynamic_cast<T*>(item.get());
        if(ptr==nullptr)
            continue;
        tree->Fill();
    }

    REQUIRE( tree->GetEntries() == expectedEntries );

    file->Write();
    file->Close();
    delete file;

    // open it again

    file = new TFile(tmpfile.filename.c_str());

    REQUIRE(file->IsOpen());

    tree = dynamic_cast<TTree*>(file->Get("tree"));
    REQUIRE((tree != nullptr));
    REQUIRE((tree->GetEntries()==expectedEntries));

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
