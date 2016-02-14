#include "base/CmdLine.h"
#include "base/Logger.h"

#include "base/WrapTFile.h"

#include "tree/TEvent.h"
#include "tree/TEventData.h"

#include "TTree.h"

using namespace ant;
using namespace std;


int main( int argc, char** argv )
{
    SetupLogger();

    TCLAP::CmdLine cmd("compare_ant_goat", ' ', "0.1");

    auto cmd_antinput = cmd.add<TCLAP::ValueArg<string>>("","ant","treeEvents from Ant",true,"","rootfile");
    auto cmd_goatinput = cmd.add<TCLAP::ValueArg<string>>("","goat","treeEvents from Goat",true,"","rootfile");

    cmd.parse(argc, argv);

    WrapTFileInput antinput(cmd_antinput->getValue());
    WrapTFileInput goatinput(cmd_goatinput->getValue());

    TTree* treeAnt;
    if(!antinput.GetObject("treeEvents", treeAnt)) {
        LOG(ERROR) << "Cannot find treeEvents in " << cmd_antinput->getValue();
        exit(EXIT_FAILURE);
    }
    LOG(INFO) << "Ant treeEvent Entries " << treeAnt->GetEntries();

    TTree* treeGoat;
    if(!goatinput.GetObject("treeEvents", treeGoat)) {
        LOG(ERROR) << "Cannot find treeEvents in " << cmd_goatinput->getValue();
        exit(EXIT_FAILURE);
    }
    LOG(INFO) << "Goat treeEvent Entries " << treeGoat->GetEntries();


    TEvent* eventAnt = new TEvent();
    TEvent* eventGoat = new TEvent();

    treeAnt->SetBranchAddress("data",addressof(eventAnt));
    treeGoat->SetBranchAddress("data",addressof(eventGoat));

    long long entryAnt = 0;
    long long entryGoat = 0;

    bool synced = false;
    long long nEvents = 0;

    while(entryAnt < treeAnt->GetEntries() && entryGoat < treeGoat->GetEntries()) {

        treeAnt->GetEntry(entryAnt);
        treeGoat->GetEntry(entryGoat);

        if(!synced) {
            const auto antEventID = eventAnt->Reconstructed->Trigger.DAQEventID;
            const auto goatEventID = eventGoat->Reconstructed->Trigger.DAQEventID;
            if(antEventID < goatEventID) {
                entryAnt++;
                continue;
            }
            else if(antEventID > goatEventID) {
                entryGoat++;
                continue;
            }
            synced = true;
            LOG(INFO) << "Synced at ant=" << entryAnt << " goat=" << entryGoat;
        }

        const TEventData& antRecon = *eventAnt->Reconstructed;
        const TEventData& goatRecon = *eventGoat->Reconstructed;

        if(antRecon.TaggerHits.size() != goatRecon.TaggerHits.size()) {
            LOG(ERROR) << "TaggerHits size mismatch at ant=" << entryAnt << " goat=" << entryGoat;
            break;
        }

        nEvents++;
        entryAnt++;
        entryGoat++;
    }

    LOG(INFO) << "Compared " << nEvents << " events";
}