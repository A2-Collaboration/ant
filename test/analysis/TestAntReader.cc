#include "catch.hpp"
#include "catch_config.h"
#include "expconfig_helpers.h"

#include "analysis/input/ant/AntReader.h"

#include "tree/TEvent.h"
#include "tree/TEventData.h"

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

void dotest_read_unpacker();

TEST_CASE("AntReader: Read from unpacker", "[analysis]") {
    test::EnsureSetup();
    dotest_read_unpacker();
}


void dotest_read_unpacker() {
    auto unpacker = Unpacker::Get(string(TEST_BLOBS_DIRECTORY)+"/Acqu_oneevent-big.dat.xz");
    auto reconstruct = std_ext::make_unique<Reconstruct>();
    AntReader reader(nullptr, move(unpacker), move(reconstruct));

    REQUIRE((reader.GetFlags() & reader_flag_t::IsSource));

    unsigned nEvents = 0;
    unsigned nCandidates = 0;
    unsigned nSlowControls = 0;
    while(true) {
        event_t event;

        if(!reader.ReadNextEvent(event))
            break;

        nEvents++;
        nCandidates += event.Reconstructed().Candidates.size();
        nSlowControls += event.Reconstructed().SlowControls.size();
    }


    REQUIRE(nEvents==221);
    REQUIRE(nSlowControls == 8);
    REQUIRE(nCandidates == 864);

}
