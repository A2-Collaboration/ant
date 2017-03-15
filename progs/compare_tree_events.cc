#include "tclap/CmdLine.h"
#include "base/Logger.h"

#include "base/WrapTFile.h"

#include "tree/TEvent.h"
#include "tree/TEventData.h"

#include "TTree.h"
#include "TRint.h"
#include "TH1D.h"
#include "TH2D.h"


#include "analysis/plot/HistogramFactory.h"
#include "analysis/plot/root_draw.h"

#include <limits>

using namespace ant;
using namespace std;
using namespace ant::analysis;
volatile bool interrupt = false;


int main( int argc, char** argv )
{
    SetupLogger();

    signal(SIGINT, [] (int) {
        cout << ">>> Interrupted" << endl;
        interrupt = true;
    });

    TCLAP::CmdLine cmd("compare_tree_events", ' ', "0.1");

    auto cmd_input1 = cmd.add<TCLAP::ValueArg<string>>("","tree1","treeEvents 1",true,"","rootfile");
    auto cmd_input2 = cmd.add<TCLAP::ValueArg<string>>("","tree2","treeEvents 2",true,"","rootfile");
    auto cmd_maxevents = cmd.add<TCLAP::MultiArg<int>>("m","maxevents","Process only max events",false,"maxevents");


    cmd.parse(argc, argv);

    WrapTFileInput input1(cmd_input1->getValue());
    WrapTFileInput input2(cmd_input2->getValue());

    TTree* tree1;
    if(!input1.GetObject("treeEvents", tree1)) {
        LOG(ERROR) << "Cannot find treeEvents in " << cmd_input1->getValue();
        exit(EXIT_FAILURE);
    }
    LOG(INFO) << "First treeEvent Entries " << tree1->GetEntries();

    TTree* tree2;
    if(!input2.GetObject("treeEvents", tree2)) {
        LOG(ERROR) << "Cannot find treeEvents in " << cmd_input2->getValue();
        exit(EXIT_FAILURE);
    }
    LOG(INFO) << "Second treeEvent Entries " << tree2->GetEntries();


    TEvent* event1 = new TEvent();
    TEvent* event2 = new TEvent();

    tree1->SetBranchAddress("data",addressof(event1));
    tree2->SetBranchAddress("data",addressof(event2));

    long long entry1 = 0;
    long long entry2 = 0;

    long long maxevents = cmd_maxevents->isSet()
            ? cmd_maxevents->getValue().back()
            :  numeric_limits<long long>::max();
    long long nEvents = 0;


    while(entry1 < tree1->GetEntries() && entry2 < tree2->GetEntries()) {

        if(interrupt)
            break;

        tree1->GetEntry(entry1);
        tree2->GetEntry(entry2);

        const TEventData& recon1 = event1->Reconstructed();
        const TEventData& recon2 = event2->Reconstructed();


        if(recon1.ID != recon2.ID) {
            LOG(INFO) << "ID mismatch ID1=" << recon1.ID << " ID2=" << recon2.ID;
            break;
        }

        if(recon1.Trigger.DAQEventID != recon2.Trigger.DAQEventID) {
            LOG(INFO) << "DAQEventID mismatch ID1=" << recon1.Trigger.DAQEventID << " ID2=" << recon2.Trigger.DAQEventID;
            break;
        }

        nEvents++;
        if(nEvents==maxevents)
            break;
        entry1++;
        entry2++;
    }

    LOG(INFO) << "Compared " << nEvents << " events";
}