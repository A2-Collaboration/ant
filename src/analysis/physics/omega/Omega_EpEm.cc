#include "Omega_EpEm.h"

#include "base/Logger.h"

#include "plot/root_draw.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;


Omega_EpEm::Omega_EpEm(const string &name, OptionsPtr opts) :
    Physics(name, opts)
{
    BinSettings bins_nClusters(20);
    BinSettings energy_binning(250);

    // so far just a copy of Tutorial class
    h_nClusters = HistFac.makeTH1D("Number of Clusters", // title
                                   "nClusters","#",      // xlabel, ylabel
                                   bins_nClusters,       // our binnings, may write directly BinSettings(10) here
                                   "h_nClusters"         // ROOT object name, auto-generated if omitted
                                   );

    h_nClusters_pr = HistFac.makeTH1D("Number of Clusters - prompt-random",
                                      "nClusters","#",
                                      bins_nClusters,
                                      "h_nClusters_pr"
                                      );

    promptrandom.AddPromptRange({ -7,   7}); // in nanoseconds
    promptrandom.AddRandomRange({-50, -10});
    promptrandom.AddRandomRange({ 10,  50});

    t.CreateBranches(HistFac.makeTTree("t"));

}

void Omega_EpEm::ProcessEvent(const TEvent& event, manager_t& manager)
{
    for(auto& taggerhit : event.Reconstructed().TaggerHits) {
        promptrandom.SetTaggerHit(taggerhit.Time);
        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;
        h_nClusters_pr->Fill(event.Reconstructed().Clusters.size(), promptrandom.FillWeight());

        t.TaggW = promptrandom.FillWeight();
        t.nClusters = event.Reconstructed().Clusters.size();
        t.Tree->Fill();
    }
    h_nClusters->Fill(event.Reconstructed().Clusters.size());

    const auto& cands = event.Reconstructed().Candidates;
    TCandidatePtrList cands_taps;
    TCandidatePtrList cands_cb;


    b_CBSumVetoE = 0;
    for(const auto& p : cands.get_iter()) {
        if(p->Detector & Detector_t::Any_t::TAPS_Apparatus) {
            cands_taps.emplace_back(p);
        }
        else if(p->Detector & Detector_t::Any_t::CB_Apparatus) {
            cands_cb.emplace_back(p);
            b_CBSumVetoE += p->VetoEnergy;
        }
    }
    b_nTAPS = cands_taps.size();

    b_nCB = cands_cb.size();


}

void Omega_EpEm::ShowResult()
{
    // ant::canvas nice wrapper around TCanvas
    ant::canvas(GetName()+": Basic plots")
            << h_nClusters
            << h_nClusters_pr
            << TTree_drawable(t.Tree, "nClusters >> (20,0,20)", "TaggW")
            << endc; // actually draws the canvas
}

AUTO_REGISTER_PHYSICS(Omega_EpEm)
