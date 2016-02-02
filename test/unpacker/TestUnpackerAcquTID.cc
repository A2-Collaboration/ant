#include "catch.hpp"
#include "catch_config.h"
#include "expconfig_helpers.h"

#include "Unpacker.h"
#include "expconfig/ExpConfig.h"

#include "tree/TEvent.h"
#include "tree/TEventData.h"

#include <iostream>
#include <string>

using namespace std;
using namespace ant;

void dotest();

TEST_CASE("Test UnpackerAcqu: TID", "[unpacker]") {
    test::EnsureSetup();
    dotest();
}

TID GetFirstID(const string& filename) {
    auto unpacker = Unpacker::Get(filename);

    auto firstevent = unpacker->NextEvent();
    if(firstevent == nullptr)
        throw runtime_error(string("Didn't get first event from ")+filename);

    return firstevent->Reconstructed->ID;
}

void dotest() {
    auto earlier = GetFirstID(string(TEST_BLOBS_DIRECTORY)+"/Acqu_headeronly-big.dat.xz");

    auto info_1 = GetFirstID(string(TEST_BLOBS_DIRECTORY)+"/Acqu_headeronly-METMEST_1.dat.xz");
    auto info_2 = GetFirstID(string(TEST_BLOBS_DIRECTORY)+"/Acqu_headeronly-METMEST_2.dat.xz");
    auto info_3 = GetFirstID(string(TEST_BLOBS_DIRECTORY)+"/Acqu_headeronly-METMEST_3.dat.xz");
    auto info_4 = GetFirstID(string(TEST_BLOBS_DIRECTORY)+"/Acqu_headeronly-METMEST_4.dat.xz");

    auto later = GetFirstID(string(TEST_BLOBS_DIRECTORY)+"/Acqu_headeronly-small.dat.xz");

    REQUIRE(earlier < info_1);
    REQUIRE(info_1 < info_2);
    REQUIRE(info_2 < info_3);
    REQUIRE(info_3 < info_4);
    REQUIRE(info_4 < later);

}
