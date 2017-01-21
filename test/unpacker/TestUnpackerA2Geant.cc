#include "catch.hpp"
#include "catch_config.h"
#include "expconfig_helpers.h"

#include "Unpacker.h"
#include "tree/TEvent.h"
#include "tree/TEventData.h"

#include <string>

using namespace std;
using namespace ant;

void dotest_readall();
void dotest_single();

TEST_CASE("UnpackerA2Geant: Read all", "[unpacker]") {
    dotest_readall();
}

TEST_CASE("UnpackerA2Geant: Single Event", "[unpacker]") {
    dotest_single();
}

struct inspect_TEvent : TEvent {
    inspect_TEvent(TEvent event) :
        TEvent(move(event)) {}
    bool HasReconstructed() const {
        return reconstructed!=nullptr;
    }
    bool HasMCTrue() const {
        return mctrue!=nullptr;
    }
};

void dotest_readall() {
    test::EnsureSetup();

    std::unique_ptr<Unpacker::Module> unpacker;
    REQUIRE_NOTHROW(unpacker = ant::Unpacker::Get(string(TEST_BLOBS_DIRECTORY)+"/Geant_with_TID.root"));
    REQUIRE(unpacker != nullptr);

    unsigned nEvents = 0;
    while(auto event = unpacker->NextEvent()) {
        nEvents++;
        REQUIRE(inspect_TEvent(move(event)).HasReconstructed());
    }

    REQUIRE(nEvents==100);
}

void dotest_single() {
    std::unique_ptr<Unpacker::Module> unpacker;
    REQUIRE_NOTHROW(unpacker = ant::Unpacker::Get(string(TEST_BLOBS_DIRECTORY)+"/Geant_with_TID.root"));
    REQUIRE(unpacker != nullptr);

    auto event = unpacker->NextEvent();
    REQUIRE(event);

    inspect_TEvent e(move(event));
    REQUIRE(e.HasMCTrue());
    REQUIRE(e.HasReconstructed());

    TEventData& recon = e.Reconstructed();
    TEventData& mctrue = e.MCTrue();

    REQUIRE(mctrue.Target.Vertex.z == Approx(-3.50488));
    REQUIRE(std::isnan(recon.Target.Vertex.z));
}