
#include "Unpacker.h"


#include "catch.hpp"
#include "test_config.h"
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
  REQUIRE_NOTHROW(ant::Unpacker::Get(filename));
}
