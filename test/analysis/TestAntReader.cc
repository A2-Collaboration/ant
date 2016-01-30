#include "catch.hpp"
#include "catch_config.h"

#include "analysis/input/ant/AntReader.h"

#include "tree/TEvent.h"

#include "expconfig/ExpConfig.h"
#include "unpacker/Unpacker.h"
#include "reconstruct/Reconstruct.h"

#include "base/WrapTFile.h"
#include "base/tmpfile_t.h"

#include "TTree.h"

#include <string>
#include <iostream>

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::input;

void dotest_read();

TEST_CASE("AntReader: Simply read", "[analysis]") {
    dotest_read();
}

void dotest_read() {
    ant::ExpConfig::Setup::ManualName = "Setup_Test";
    auto unpacker = Unpacker::Get(string(TEST_BLOBS_DIRECTORY)+"/Acqu_oneevent-big.dat.xz");
    auto reconstruct = std_ext::make_unique<Reconstruct>();
    AntReader reader(nullptr, move(unpacker), move(reconstruct));

    REQUIRE(reader.IsSource());

    unsigned nEvents = 0;
    unsigned nCandidates = 0;
    unsigned nSlowControls = 0;
    while(true) {
        TEvent event;

        if(!reader.ReadNextEvent(event))
            break;

        nEvents++;
        nCandidates += event.Reconstructed->Candidates.size();
        nSlowControls += event.Reconstructed->SlowControls.size();
    }


    REQUIRE(nEvents==221);
    REQUIRE(nSlowControls == 8);
    REQUIRE(nCandidates == 822);

}


