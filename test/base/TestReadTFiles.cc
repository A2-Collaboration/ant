#include "catch.hpp"
#include "catch_config.h"

#include "base/ReadTFiles.h"

#include "unpacker/Unpacker.h"

#include "tree/THeaderInfo.h"

#include "base/WrapTFile.h"
#include "base/tmpfile_t.h"


#include "TTree.h"

#include <iostream>

using namespace std;

void dotest();
void generateInputFile(const std::string& filename);


TEST_CASE("AntReader", "[analysis]") {
    dotest();
}

void dotest() {

    ant::tmpfile_t tmp;
    ant::ReadTFiles m;

    generateInputFile(tmp.filename);

    m.OpenFile(tmp.filename);

    TTree* tree = nullptr;

    REQUIRE( m.GetObject("treeHeaderInfo", tree));
    REQUIRE(tree != nullptr);
    cout << "Found tree " << tree->GetName() << endl;
    REQUIRE(tree->GetEntries() == 1);

    m.CloseAll();
}


void generateInputFile(const string& filename) {

    auto unpacker = ant::Unpacker::Get(string(TEST_BLOBS_DIRECTORY)+"/Acqu_oneevent-small.dat.xz");

    // write some stuff to a ROOT tree
    ant::WrapTFile file(filename);
    TTree* treeHeaderInfo = file.CreateInside<TTree>("treeHeaderInfo", "treeHeaderInfo");
    ant::THeaderInfo* HeaderInfo = file.CreateInside<ant::THeaderInfo>();

    treeHeaderInfo->Branch("THeaderInfo", "ant::THeaderInfo", &HeaderInfo);

    while(auto item = unpacker->NextItem()) {
        HeaderInfo = dynamic_cast<ant::THeaderInfo*>(item.get());
        if(HeaderInfo != nullptr) {
            treeHeaderInfo->Fill();
            break;
        }
    }
}
