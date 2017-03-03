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
    gPad->SetLogz();
    h_cb = new TH2CB("h_cb","CB");
    h_cb->Draw("colz");

    cd(2);
    gPad->SetLogz();
    h_taps = new TH2TAPS("h_taps","TAPS");
    h_taps->Draw("colz");

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

    h_cb->ResetElements();
    h_taps->ResetElements();
    h_cb->ClearMarkers();
    h_taps->ClearMarkers();

    cout << ">>>> EVENT " << recon.ID << endl;
    cout << "Candidates n=" << recon.Candidates.size() << endl;
    cout << "Clusters n=" << recon.Clusters.size() << endl;
    for(auto& cand : recon.Candidates) {
        cout << cand << endl;
        // show only clusters belonging to candidates (above cluster threshold)
        auto caloCluster = cand.FindCaloCluster();
        auto h = caloCluster->DetectorType == Detector_t::Type_t::CB ?
                     dynamic_cast<TH2Crystals*>(h_cb) : h_taps;
        h->CreateMarker(caloCluster->CentralElement);
        cout << *caloCluster << endl;
    }
    cout << endl;

    // show all energies which was clustered
    for(auto& cluster : recon.Clusters) {
        auto h = cluster.DetectorType == Detector_t::Type_t::CB ?
                     dynamic_cast<TH2Crystals*>(h_cb) : h_taps;
        for(auto& hit : cluster.Hits) {
            h->SetElement(hit.Channel, hit.Energy);
        }
    }

    cd(1);
    gPad->Modified();
    gPad->Update();

    cd(2);
    gPad->Modified();
    gPad->Update();
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
