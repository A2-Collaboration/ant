#include "catch.hpp"
#include "catch_config.h"

#include "Unpacker.h"

#include "expconfig_helpers.h"

#include <iostream>
#include <string>

using namespace std;

void dotest(const string& filename);

TEST_CASE("Test UnpackerAcqu: headeronly, big record length", "[unpacker]") {
    dotest(string(TEST_BLOBS_DIRECTORY)+"/Acqu_headeronly-big.dat.xz");
}

TEST_CASE("Test UnpackerAcqu: one event block, big record length", "[unpacker]") {
    dotest(string(TEST_BLOBS_DIRECTORY)+"/Acqu_oneevent-big.dat.xz");
}

TEST_CASE("Test UnpackerAcqu: headeronly, small record length", "[unpacker]") {
    dotest(string(TEST_BLOBS_DIRECTORY)+"/Acqu_headeronly-small.dat.xz");
}

TEST_CASE("Test UnpackerAcqu: one event block, small record length", "[unpacker]") {
    dotest(string(TEST_BLOBS_DIRECTORY)+"/Acqu_oneevent-small.dat.xz");
}

TEST_CASE("Test UnpackerAcqu: Mk1, three blocks", "[unpacker]") {
    dotest(string(TEST_BLOBS_DIRECTORY)+"/AcquMk1_problematic.dat.gz");
}

TEST_CASE("TestUnpackerAcqu: Mk2, big buffer size of 32*0x8000", "[unpacker]") {
	dotest(string(TEST_BLOBS_DIRECTORY)+"/Acqu_headeronly_increasedBufferSize.dat.xz");
}

void dotest(const string &filename) {
    ant::test::EnsureSetup();
    // this simply tries to open the file
    REQUIRE_NOTHROW(ant::Unpacker::Get(filename));
}
