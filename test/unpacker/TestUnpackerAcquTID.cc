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

TEST_CASE("Test UnpackerAcqu: TID", "[unpacker]") {
    dotest();
}

THeaderInfo GetHeaderInfo(const string& filename) {
    ExpConfig::Setup::ManualName = "Setup_Test";
    auto unpacker = ant::Unpacker::Get(filename);

    while(auto item = unpacker->NextItem()) {
        auto HeaderInfo = dynamic_cast<THeaderInfo*>(item.get());
        if(HeaderInfo != nullptr) {
            return *HeaderInfo;
        }
    }
    throw runtime_error(string("Cannot find header info in ")+filename);
}

void dotest() {
    THeaderInfo earlier = GetHeaderInfo(string(TEST_BLOBS_DIRECTORY)+"/Acqu_headeronly-big.dat.xz");

    THeaderInfo info_1 = GetHeaderInfo(string(TEST_BLOBS_DIRECTORY)+"/Acqu_headeronly-METMEST_1.dat.xz");
    THeaderInfo info_2 = GetHeaderInfo(string(TEST_BLOBS_DIRECTORY)+"/Acqu_headeronly-METMEST_2.dat.xz");
    THeaderInfo info_3 = GetHeaderInfo(string(TEST_BLOBS_DIRECTORY)+"/Acqu_headeronly-METMEST_3.dat.xz");
    THeaderInfo info_4 = GetHeaderInfo(string(TEST_BLOBS_DIRECTORY)+"/Acqu_headeronly-METMEST_4.dat.xz");

    THeaderInfo later = GetHeaderInfo(string(TEST_BLOBS_DIRECTORY)+"/Acqu_headeronly-small.dat.xz");

    cout << earlier << endl;
    cout << info_1 << endl;
    cout << info_2 << endl;
    cout << info_3 << endl;
    cout << info_4 << endl;
    cout << later << endl;

    REQUIRE(earlier.ID < info_1.ID);
    REQUIRE(info_1.ID < info_2.ID);
    REQUIRE(info_2.ID < info_3.ID);
    REQUIRE(info_3.ID < info_4.ID);
    REQUIRE(info_4.ID < later.ID);

}
