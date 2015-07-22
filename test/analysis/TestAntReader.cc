#include "catch.hpp"
#include "catch_config.h"

#include "analysis/input/ant/AntReader.h"
#include "analysis/data/Event.h"

#include "tree/THeaderInfo.h"
#include "tree/TDetectorRead.h"
#include "tree/TEvent.h"

#include "expconfig/ExpConfig.h"
#include "unpacker/Unpacker.h"
#include "reconstruct/Reconstruct.h"

#include "base/TFileWrapper.h"
#include "base/tmpfile_t.h"

#include "TTree.h"

#include <string>
#include <iostream>

using namespace std;
using namespace ant;
using namespace ant::input;

void dotest();
void generateInputFile(const std::string& filename);

TEST_CASE("AntReader", "[analysis]") {
    dotest();
}


void dotest() {
    tmpfile_t tmp;
    generateInputFile(tmp.filename);

    ant::input::AntReader reader;

    reader.AddInputFile(tmp.filename);

    reader.Initialize();

    unsigned int nEvents = 0;
    while(reader.hasData()) {
        auto event = reader.ReadNextEvent();
        REQUIRE(event != nullptr);
        /// \bug Fix missing tracks
        //REQUIRE(event->Reconstructed().Tracks().size()>0);
        nEvents++;
    }
    REQUIRE(nEvents==222);
}


void generateInputFile(const string& filename) {

    auto unpacker = Unpacker::Get(string(TEST_BLOBS_DIRECTORY)+"/Acqu_oneevent-big.dat.xz");

    // write some stuff to a ROOT tree
    TFileWrapper file(filename);

    TTree* treeEvent = new TTree("treeEvent", "treeEvent");
    TEvent* Event = new TEvent();
    treeEvent->Branch("Event", "ant::TEvent", &Event);


    unique_ptr<Reconstruct> reconstruct;

    while(auto item = unpacker->NextItem()) {

        auto HeaderInfo = dynamic_cast<THeaderInfo*>(item.get());
        if(HeaderInfo != nullptr) {
            reconstruct = std_ext::make_unique<Reconstruct>(*HeaderInfo);
            continue;
        }

        auto DetectorRead = dynamic_cast<TDetectorRead*>(item.get());
        if(DetectorRead != nullptr) {

            if(reconstruct) {
                auto event_ptr = reconstruct->DoReconstruct(*DetectorRead);
                Event = event_ptr.get();
                treeEvent->Fill();
            }
        }
    }
}
