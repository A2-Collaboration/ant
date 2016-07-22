#include "catch.hpp"
#include "catch_config.h"

#include "base/WrapTFile.h"
#include "tree/TEvent.h"
#include "tree/TEventData.h"

#include "analysis/input/pluto/PlutoReader.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::input;

void dotest();

TEST_CASE("Fitter: TODO", "[analysis]") {
    dotest();
}


void dotest() {
    auto rootfile = make_shared<WrapTFileInput>(string(TEST_BLOBS_DIRECTORY)+"/Pluto_Etap2g.root");
    PlutoReader reader(rootfile);

    REQUIRE_FALSE(reader.IsSource());

    unsigned nEvents = 0;
    while(true) {
        TEvent event;

        if(!reader.ReadNextEvent(event))
            break;

        nEvents++;
    }

    REQUIRE(nEvents==500);
}
