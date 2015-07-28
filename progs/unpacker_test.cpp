
#include <iostream>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include "unpacker/Unpacker.h"
#include "expconfig/ExpConfig.h"

#include "tree/THeaderInfo.h"
#include "tree/TUnpackerMessage.h"
#include "tree/TSlowControl.h"
#include "tree/TDetectorRead.h"
#include "tree/TEvent.h"

#include "unpacker/RawFileReader.h"

#include "reconstruct/Reconstruct.h"

#include "base/Logger.h"
#include "base/Format.h"
#include "base/WrapTFile.h"
#include "base/CmdLine.h"

#include <chrono>

#include "TTree.h"

using namespace std;
using namespace ant;

bool running = true;

void myCrashHandler(int sig) {
    if(sig == SIGINT) {
        running = false;
        return;
    }
    // FOLLOWING LINE IS ABSOLUTELY NEEDED AT THE END IN ORDER TO ABORT APPLICATION
    el::Helpers::crashAbort(sig);
}

int main(int argc, char* argv[]) {
    SetupLogger();

    TCLAP::CmdLine cmd("unpacker_test", ' ', "0.01");
    auto cmd_input  = cmd.add<TCLAP::ValueArg<string>>("i","input","Unpacker input file (Acqu/A2Geant)",false,"scratch/steffen-pi0-cluster-geant.root","string");
    auto cmd_setup  = cmd.add<TCLAP::ValueArg<string>>("","setup","Manually choose setup",false,"Setup_2014_EtaPrime","string");
    auto cmd_output = cmd.add<TCLAP::ValueArg<string>>("o","output","Output file",false,"unpacker_test_out.root","string");
    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"int");
    auto cmd_noreads = cmd.add<TCLAP::SwitchArg>("","noreads","Do not write detector reads to output",false);
    auto cmd_noevents = cmd.add<TCLAP::SwitchArg>("","noevents","Do not write events to output",false);

    cmd.parse(argc, argv);
    if(cmd_verbose->isSet()) {
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());
    }
    el::Helpers::setCrashHandler(myCrashHandler);


    ExpConfig::ManualSetupName = cmd_setup->getValue(); // if non-empty, disables automatic setup determination
    RawFileReader::OutputPerformanceStats = 5;

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    //auto unpacker = Unpacker::Get("scratch/CBTaggTAPS_9227.dat");
    //auto unpacker = Unpacker::Get("scratch/CBTaggTAPS_7892.dat");
    //auto unpacker = Unpacker::Get("scratch/CBTaggTAPS_5711.dat.xz");
    //auto unpacker = Unpacker::Get("scratch/oneevent-small.dat");
    auto unpacker = Unpacker::Get(cmd_input->getValue());
    unsigned nReads = 0;

    // write some stuff to a ROOT tree
    WrapTFile file(cmd_output->getValue());

    TTree* treeHeaderInfo = new TTree("treeHeaderInfo", "treeHeaderInfo");
    THeaderInfo* HeaderInfo = new THeaderInfo();
    treeHeaderInfo->Branch("THeaderInfo", "ant::THeaderInfo", &HeaderInfo);

    TTree* treeUnpackerMessage = new TTree("treeUnpackerMessage", "treeUnpackerMessage");
    TUnpackerMessage* UnpackerMessage = new TUnpackerMessage();
    treeUnpackerMessage->Branch("UnpackerMessage", "ant::TUnpackerMessage", &UnpackerMessage);

    TTree* treeSlowControl = new TTree("treeSlowControl", "treeSlowControl");
    TSlowControl* SlowControl = new TSlowControl();
    treeSlowControl->Branch("SlowControl", "ant::TSlowControl", &SlowControl);

    TTree* treeDetectorRead = new TTree("treeDetectorRead", "treeDetectorRead");
    TDetectorRead* DetectorRead = new TDetectorRead();
    treeDetectorRead->Branch("DetectorRead", "ant::TDetectorRead", &DetectorRead);

    TTree* treeEvent = new TTree("treeEvent", "treeEvent");
    TEvent* Event = new TEvent();
    treeEvent->Branch("Event", "ant::TEvent", &Event);


    unique_ptr<Reconstruct> reconstruct;

    while(auto item = unpacker->NextItem()) {
        if(!running)
            break;

        //cout << *item << endl;
        HeaderInfo = dynamic_cast<THeaderInfo*>(item.get());
        if(HeaderInfo != nullptr) {
            treeHeaderInfo->Fill();
            reconstruct = std_ext::make_unique<Reconstruct>(*HeaderInfo);
            continue;
        }
        UnpackerMessage = dynamic_cast<TUnpackerMessage*>(item.get());
        if(UnpackerMessage != nullptr) {
            treeUnpackerMessage->Fill();
            continue;
        }
        SlowControl = dynamic_cast<TSlowControl*>(item.get());
        if(SlowControl != nullptr) {
            treeSlowControl->Fill();
            continue;
        }
        DetectorRead = dynamic_cast<TDetectorRead*>(item.get());
        if(DetectorRead != nullptr) {

            if(reconstruct) {
                auto event_ptr = reconstruct->DoReconstruct(*DetectorRead);
                Event = event_ptr.get();
                if(!cmd_noevents->isSet())
                    treeEvent->Fill();
            }

            if(!cmd_noreads->isSet())
                treeDetectorRead->Fill();

            nReads++;
            if(nReads % 3000 == 0) {
                LOG(INFO) << "Unpacked " << nReads << " detector reads";
            }
            continue;
        }
    }

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    cout << "Analyzed " << nReads << " reads, speed " << nReads/elapsed_seconds.count() << " Reads/s" << endl;

    return EXIT_SUCCESS;
}
