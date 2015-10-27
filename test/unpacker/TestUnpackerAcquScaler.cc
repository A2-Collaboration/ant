#include "catch.hpp"
#include "catch_config.h"

#include "Unpacker.h"
#include "expconfig/ExpConfig.h"

#include "tree/THeaderInfo.h"
#include "tree/TDetectorRead.h"
#include "tree/TSlowControl.h"

#include <iostream>
#include <string>

using namespace std;
using namespace ant;

void dotest();

TEST_CASE("Test UnpackerAcqu: Scaler block", "[unpacker]") {
    dotest();
}

void dotest() {
    ant::ExpConfig::Setup::ManualName = "Setup_Test";
    auto unpacker = ant::Unpacker::Get(string(TEST_BLOBS_DIRECTORY)+"/Acqu_scalerblock.dat.xz");

    unsigned nSlowControls = 0;
    unsigned nReads = 0;
    unsigned nHits = 0;

    bool taggerScalerBlockFound = false;

    while(auto item = unpacker->NextItem()) {
        auto DetectorRead = dynamic_cast<TDetectorRead*>(item.get());
        if(DetectorRead != nullptr) {
            nReads++;
            nHits += DetectorRead->Hits.size();
        }
        auto SlowControl = dynamic_cast<TSlowControl*>(item.get());
        if(SlowControl != nullptr) {
            nSlowControls++;
            if(SlowControl->Name == "TaggerScalers") {
                taggerScalerBlockFound = true;
                // we know the file is extracted from an EPT run...
                REQUIRE(SlowControl->Payload_Int.size() == 47);
            }
        }
    }

    REQUIRE(nSlowControls == 15);
    REQUIRE(nReads == 211);
    REQUIRE(nHits == 28456);
    REQUIRE(taggerScalerBlockFound);
}
