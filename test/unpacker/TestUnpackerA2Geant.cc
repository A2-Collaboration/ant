#include "catch.hpp"
#include "catch_config.h"

#include "Unpacker.h"
#include "expconfig/ExpConfig.h"

#include <string>

using namespace std;

TEST_CASE("Test UnpackerA2Geant", "[unpacker]") {
    // try to open the file, but we need some setup name for this...
    ant::ExpConfig::ManualSetupName = "Setup_2014_07_EPT_Prod";
    REQUIRE_NOTHROW(ant::Unpacker::Get(string(TEST_BLOBS_DIRECTORY)+"/Geant_with_TID.root"));
}

