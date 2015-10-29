#include "catch.hpp"
#include "catch_config.h"

#include "analysis/input/ant/AntReader.h"
#include "analysis/data/Event.h"

#include "tree/THeaderInfo.h"
#include "tree/TDetectorRead.h"
#include "tree/TEvent.h"
#include "tree/UnpackerReader.h"

#include "expconfig/ExpConfig.h"
#include "unpacker/Unpacker.h"
#include "reconstruct/Reconstruct.h"

#include "base/WrapTFile.h"
#include "base/tmpfile_t.h"

#include "TTree.h"

#include <string>
#include <iostream>

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::input;

void dotest_read();
void dotest_writeread(bool no_reconstruct, bool uncalibratedDetectorReads, bool calibratedDetectorReads);

TEST_CASE("AntReader: Simply read", "[analysis]") {
    dotest_read();
}
TEST_CASE("AntReader: Write and read", "[analysis]") {
    dotest_writeread(false, false, false);
}
TEST_CASE("AntReader: Write and read with uncalibrated detectorreads", "[analysis]") {
    dotest_writeread(false, true, false);
}
TEST_CASE("AntReader: Write and read with calibrated detectorreads", "[analysis]") {
    dotest_writeread(false, false, true);
}
TEST_CASE("AntReader: Write and read (no reconstruct)", "[analysis]") {
    dotest_writeread(true, false, false);
}
TEST_CASE("AntReader: Write and read with uncalibrated detectorreads (no reconstruct)", "[analysis]") {
    dotest_writeread(true, true, false);
}
TEST_CASE("AntReader: Write and read with calibrated detectorreads (no reconstruct)", "[analysis]") {
    dotest_writeread(true, false, true);
}

void drive_reader(unique_ptr<AntReader> unpacker_reader, bool no_events_expected = false) {
    unsigned nEvents = 0;
    unsigned nCandidates = 0;
    unsigned nSlowControls = 0;
    while(true) {
        bool flag = false;
        data::Event event;

        if(unpacker_reader->ReadNextEvent(event)) {
            nEvents++;
            nCandidates += event.Reconstructed().Candidates().size();
            flag = true;
        }

        auto slowcontrol = unpacker_reader->ReadNextSlowControl();
        if(slowcontrol != nullptr) {
            nSlowControls++;
            flag = true;
        }

        if(!flag)
            break;
    }

    if(no_events_expected) {
        REQUIRE(nEvents==0);
    }
    else {
        REQUIRE(nEvents==221);
        REQUIRE(nSlowControls == 8);
        REQUIRE(nCandidates == 822);
    }
}

void dotest_read() {
    ant::ExpConfig::Setup::ManualName = "Setup_Test";
    auto unpacker = Unpacker::Get(string(TEST_BLOBS_DIRECTORY)+"/Acqu_oneevent-big.dat.xz");
    auto reconstruct = std_ext::make_unique<Reconstruct>();
    auto unpacker_reader = std_ext::make_unique<AntReader>(move(unpacker), move(reconstruct));

    REQUIRE(unpacker_reader->IsSource());

    drive_reader(move(unpacker_reader));
}

void dotest_writeread(
        bool no_reconstruct,
        bool uncalibratedDetectorReads,
        bool calibratedDetectorReads) {

    // FIRST STAGE (OUTPUT)
    ant::ExpConfig::Setup::ManualName = "Setup_Test";
    auto unpacker = Unpacker::Get(string(TEST_BLOBS_DIRECTORY)+"/Acqu_oneevent-big.dat.xz");
    auto reconstruct = std_ext::make_unique<Reconstruct>();
    auto unpacker_reader = std_ext::make_unique<AntReader>(move(unpacker), move(reconstruct));

    tmpfile_t output;
    REQUIRE_THROWS_AS(unpacker_reader->EnableUnpackerWriter(output.filename, true, true), DataReader::Exception);
    REQUIRE_NOTHROW(unpacker_reader->EnableUnpackerWriter(output.filename, uncalibratedDetectorReads, calibratedDetectorReads));

    // drive the reader, generates the outputfile
    drive_reader(move(unpacker_reader));

    // SECOND STAGE (READ BACK)

    unpacker_reader = nullptr;
    auto rootfiles = make_shared<WrapTFileInput>();
    REQUIRE_NOTHROW(rootfiles->OpenFile(output.filename));
    auto unpackerFile = std_ext::make_unique<tree::UnpackerReader>(rootfiles);

    REQUIRE(unpackerFile->OpenInput());

    // construct another unpacker reader from the ROOT file
    reconstruct = no_reconstruct ? nullptr : std_ext::make_unique<Reconstruct>();
    unpacker_reader = std_ext::make_unique<AntReader>(move(unpackerFile), move(reconstruct));

    // drive the reader again
    // if we ran with no reconstruct in the second stage, we expect no events
    // when reading detectorreads
    bool no_events_expected = no_reconstruct && (uncalibratedDetectorReads || calibratedDetectorReads);
    drive_reader(move(unpacker_reader), no_events_expected);
}
