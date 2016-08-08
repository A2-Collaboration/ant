#include "DisplayClustering.h"



#include "expconfig/ExpConfig.h"
#include "reconstruct/Reconstruct.h"
#include "tree/TAntHeader.h"
#include "tree/TEvent.h"
#include "tree/TEventData.h"
#include "base/Logger.h"

#include "root-addons/cbtaps_display/TH2CB.h"
#include "root-addons/cbtaps_display/TH2TAPS.h"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

#include <memory>

using namespace ant;
using namespace std;

namespace ant {
namespace detail {

struct Implementation {
    std::unique_ptr<Reconstruct_traits> reconstruct;
};

}} // namespace ant::detail

DisplayClustering::DisplayClustering(TFile* file)
{
    SetupLogger();

    if(!file) {
        auto files = gROOT->GetListOfFiles();
        if(files->GetEntries() != 1) {
            LOG(ERROR) << "Can only work on 1 file";
            return;
        }
        file = dynamic_cast<TFile*>(files->At(0));
    }

    TAntHeader* AntHeader = nullptr;
    file->GetObject("AntHeader", AntHeader);
    if(!AntHeader) {
        LOG(ERROR) << "File did not contain AntHeader";
        return;
    }

    auto setup = ExpConfig::Setup::Get(AntHeader->SetupName);
    if(!setup) {
        LOG(ERROR) << "Setup not found";
        return;
    }

    file->GetObject("treeEvents", treeEvents);
    if(!treeEvents) {
        LOG(ERROR) << "Could not find treeEvents in file";
        return;
    }

    const auto res = treeEvents->SetBranchAddress("data", addressof(eventPtr));
    if(res != TTree::kMatch) {
        LOG(ERROR) << "Could not access branch 'data' in treeEvents";
        return;
    }

    Divide(2,1);

    cd(1);
    h_cb = new TH2CB("h_cb","CB");
    h_cb->Draw();

    cd(2);
    h_taps = new TH2TAPS("h_taps","TAPS");
    h_taps->Draw();

    impl = new detail::Implementation();

    impl->reconstruct = std_ext::make_unique<Reconstruct>();

    Display();
}

DisplayClustering::~DisplayClustering()
{
    if(impl)
        delete impl;
}

void DisplayClustering::Display()
{
    if(!impl)
        return;
    if(curr_entry >= treeEvents->GetEntries())
        return;
    treeEvents->GetEntry(curr_entry);

    auto& recon = eventPtr->Reconstructed();
    impl->reconstruct->DoReconstruct(recon);

    cout << ">>>> EVENT " << recon.ID << endl;
    cout << "Candidates n=" << recon.Candidates.size() << endl;
    for(auto& cand : recon.Candidates) {
        cout << cand << endl;
    }
    cout << endl;


}

void DisplayClustering::Next()
{
    curr_entry++;
    if(curr_entry >= treeEvents->GetEntries())
        curr_entry = treeEvents->GetEntries()-1;
    Display();
}

void DisplayClustering::Prev()
{
    curr_entry--;
    if(curr_entry < 0)
        curr_entry = 0;
    Display();
}

void DisplayClustering::HandleInput(EEventType button, Int_t x, Int_t y)
{
    if(button == kKeyPress) {
        if(char(x) == 'd') {
            Next();
        } else if( char(x) == 'a') {
            Prev();
        }
    }
    TCanvas::HandleInput(button, x ,y);
}
