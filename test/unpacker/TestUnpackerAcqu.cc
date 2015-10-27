#include "catch.hpp"
#include "catch_config.h"

#include "Unpacker.h"
#include "expconfig/ExpConfig.h"

#include <iostream>
#include <string>

using namespace std;

void dotest(const string& filename);

TEST_CASE("Test UnpackerAcqu: headeronly, big record length", "[unpacker]") {
    dotest(string(TEST_BLOBS_DIRECTORY)+"/Acqu_headeronly-big.dat.xz");
}

TEST_CASE("Test UnpackerAcqu: one event, big record length", "[unpacker]") {
    dotest(string(TEST_BLOBS_DIRECTORY)+"/Acqu_oneevent-big.dat.xz");
}

TEST_CASE("Test UnpackerAcqu: headeronly, small record length", "[unpacker]") {
    dotest(string(TEST_BLOBS_DIRECTORY)+"/Acqu_headeronly-small.dat.xz");
}

TEST_CASE("Test UnpackerAcqu: one event, small record length", "[unpacker]") {
    dotest(string(TEST_BLOBS_DIRECTORY)+"/Acqu_oneevent-small.dat.xz");
}


void dotest(const string &filename) {
    ant::ExpConfig::Setup::ManualName = "Setup_Test";
    // this simply tries to open the file
    REQUIRE_NOTHROW(ant::Unpacker::Get(filename));
}
