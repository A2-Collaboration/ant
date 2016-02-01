#include "catch.hpp"
#include "catch_config.h"
#include "expconfig_helpers.h"

#include "Unpacker.h"
#include "tree/TEvent.h"

#include <string>

using namespace std;
using namespace ant;

void dotest();

TEST_CASE("Test UnpackerA2Geant", "[unpacker]") {
    dotest();
}

void dotest() {
    test::EnsureSetup();

    std::unique_ptr<Unpacker::Module> unpacker;
    REQUIRE_NOTHROW(unpacker = ant::Unpacker::Get(string(TEST_BLOBS_DIRECTORY)+"/Geant_with_TID.root"));
    REQUIRE(unpacker != nullptr);

    unsigned nEvents = 0;
    while(auto event = unpacker->NextEvent()) {
        nEvents++;
        REQUIRE(event->Reconstructed != nullptr);
    }

    REQUIRE(nEvents==10);
}

