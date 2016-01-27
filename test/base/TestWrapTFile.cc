#include "catch.hpp"
#include "catch_config.h"

#include "base/WrapTFile.h"
#include "base/tmpfile_t.h"
#include "base/std_ext/memory.h"


#include "unpacker/Unpacker.h"
#include "expconfig/ExpConfig.h"

#include "TH1D.h"
#include "TTree.h"

using namespace std;
using namespace ant;

void dotest_rw();
void dotest_r();

TEST_CASE("WrapTFileInput", "[base]") {
    dotest_r();
}

TEST_CASE("WrapTFileOutput", "[base]") {
    dotest_rw();
}

void dotest_r() {
    /// \todo write some WrapTFileInput tests...
    /// create ROOT file manually, then inspect it with WrapTFileInput
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


