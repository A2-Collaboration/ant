#include "catch.hpp"
#include "catch_config.h"
#include "expconfig_helpers.h"

#include "Unpacker.h"

#include "ExpConfig.h"
#include "expconfig/detectors/Trigger.h"
#include "reconstruct/Reconstruct.h"

#include "tree/TEvent.h"
#include "tree/TEventData.h"

#include <iostream>
#include <string>

using namespace std;
using namespace ant;

void dotest();

TEST_CASE("Trigger Patterns", "[expconfig]") {
    dotest();
}

void dotest() {
    ant::test::EnsureSetup();

    auto trigger = ExpConfig::Setup::GetDetector<expconfig::detector::Trigger_2014>();
    REQUIRE(trigger);

    auto unpacker = Unpacker::Get(string(TEST_BLOBS_DIRECTORY)+"/Acqu_scalerblock.dat.xz");

    Reconstruct reconstruct;

    std::map<unsigned long, unsigned> possible_L2;
    while(auto event = unpacker->NextEvent()) {

        reconstruct.DoReconstruct(event.Reconstructed());

        CHECK(trigger->GetL1Pattern().to_ulong() == 0x1);
        unsigned lower_L2 = trigger->GetL2Pattern().to_ulong() & 0xff;
        CHECK(lower_L2 == trigger->GetL1Pattern().to_ulong());
        possible_L2[trigger->GetL2Pattern().to_ulong()]++;
        CHECK(trigger->GetMultiplicityPattern().count() == 0);
        CHECK(trigger->GetMultiplicityValue() == 0);
        CHECK(trigger->GetHelicityPattern().to_ulong() == 0xa208);
        CHECK(trigger->GetTriggerFiredPattern().to_ulong() == 0x1);
    }

    CHECK(possible_L2.size() == 16);
}
