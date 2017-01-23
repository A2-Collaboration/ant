#include "catch.hpp"
#include "catch_config.h"

#include "tree/TEventData.h"

#include "base/WrapTFile.h"
#include "analysis/input/event_t.h"
#include "analysis/input/pluto/PlutoReader.h"
#include "analysis/utils/MCWeighting.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;

void dotest();

TEST_CASE("Analysis: Raw Input", "[analysis]") {
    dotest();
}

void dotest() {
    input::PlutoReader reader(make_shared<WrapTFileInput>(string(TEST_BLOBS_DIRECTORY)+"/Pluto_EtapOmegaG.root"));

    utils::MCWeighting mcWeighting;

    while(true) {
        input::event_t event;
        if(!reader.ReadNextEvent(event))
            break;
        double w = mcWeighting.GetWeight(event.MCTrue().ParticleTree);
        /// \todo Improve tests...
        REQUIRE(w>=0);
    }


}
