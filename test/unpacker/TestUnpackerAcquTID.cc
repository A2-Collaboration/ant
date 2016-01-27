#include "catch.hpp"
#include "catch_config.h"

#include "Unpacker.h"
#include "expconfig/ExpConfig.h"

#include "tree/TDataRecord.h"

#include <iostream>
#include <string>

using namespace std;
using namespace ant;

void dotest();

TEST_CASE("Test UnpackerAcqu: TID", "[unpacker]") {
    dotest();
}

TID GetFirstID(const string& filename) {
    ExpConfig::Setup::ManualName = "Setup_Test";
    auto unpacker = ant::Unpacker::Get(filename);

    auto firstitem = unpacker->NextItem();
    if(firstitem == nullptr)
        throw runtime_error(string("Didn't get first item from ")+filename);

    return firstitem->ID;
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
