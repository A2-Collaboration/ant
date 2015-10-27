#include "catch.hpp"
#include "catch_config.h"

#include "base/WrapTFile.h"
#include "base/tmpfile_t.h"
#include "base/std_ext/memory.h"

#include "tree/THeaderInfo.h"

#include "unpacker/Unpacker.h"
#include "expconfig/ExpConfig.h"

#include "TH1D.h"
#include "TTree.h"

using namespace std;
using namespace ant;

void dotest_rw();
void dotest_r();
void generateInputFile(const std::string& filename);


TEST_CASE("WrapTFileInput", "[base]") {
    dotest_r();
}

TEST_CASE("WrapTFileOutput", "[base]") {
    dotest_rw();
}

void dotest_r() {
    ant::tmpfile_t tmp;
    ant::WrapTFileInput m;

    generateInputFile(tmp.filename);

    m.OpenFile(tmp.filename);

    TTree* tree = nullptr;

    REQUIRE( m.GetObject("treeHeaderInfo", tree));
    REQUIRE(tree != nullptr);
    REQUIRE(tree->GetEntries() == 1);
}




void dotest_rw() {
    tmpfile_t tmp;

    auto outfile = std_ext::make_unique<WrapTFileOutput>(tmp.filename);
    auto h2 = outfile->CreateInside<TH1D>("b","B",10,0,10);
    auto h3 = outfile->CreateInside<TH1D>("c","C",10,0,10);
    h2->Fill(3);
    h3->Fill(2);
    outfile = nullptr;

    WrapTFileInput infile(tmp.filename);
    REQUIRE(infile.GetListOf<TH1D>().size() == 2);
}

void generateInputFile(const string& filename) {
    ant::ExpConfig::Setup::ManualName = "Setup_Test";
    auto unpacker = ant::Unpacker::Get(string(TEST_BLOBS_DIRECTORY)+"/Acqu_oneevent-small.dat.xz");

    // write some stuff to a ROOT tree
    ant::WrapTFileOutput file(filename);
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

