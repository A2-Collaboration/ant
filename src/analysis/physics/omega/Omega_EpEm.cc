#include "Omega_EpEm.h"
#include "utils/Combinatorics.h"
#include "utils/ParticleTools.h"
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
    BinSettings bins_nParticles(10);
    BinSettings energy_binning(250);
    BinSettings im_binning(200);

    // so far just a copy of Tutorial class
    h_nClusters = HistFac.makeTH1D("Number of Clusters", // title
                                   "nClusters","#",      // xlabel, ylabel
                                   bins_nClusters,       // binning
                                   "h_nClusters"         // ROOT object name
                                   );

    h_nClusters_pr = HistFac.makeTH1D("Number of Clusters - prompt-random",
                                      "nClusters","#",
                                      bins_nClusters,
                                      "h_nClusters_pr"
                                      );

    h_nPhotons = HistFac.makeTH1D("Number of photons",
                                  "N","#",
                                  bins_nParticles,
                                  "N_g"
                                  );
    h_IM_2g = HistFac.makeTH1D("2 #gamma IM",
                                      "M_{#gamma #gamma} [MeV]","#",
                                      im_binning,
                                      "IM_2g"
                                      );

    // prompt random window
    promptrandom.AddPromptRange({ -7,   7}); // in ns
    promptrandom.AddRandomRange({-50, -10});
    promptrandom.AddRandomRange({ 10,  50});

    // some tree
    t.CreateBranches(HistFac.makeTTree("t"));

}

void Omega_EpEm::ProcessEvent(const TEvent& event, manager_t&)
{
    triggersimu.ProcessEvent(event);
    for(auto& taggerhit : event.Reconstructed().TaggerHits) {
        promptrandom.SetTaggerTime(triggersimu.GetCorrectedTaggerTime(taggerhit));
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

    auto recon_particles = utils::ParticleTypeList::Make(event.Reconstructed().Candidates);
    const auto& photons = recon_particles.Get(ParticleTypeDatabase::Photon);
    h_nPhotons->Fill(photons.size());
    //utils::ParticleTools::FillIMCombinations([this] (double x) {h_IM_2g->Fill(x);},  2, photons);
    utils::ParticleTools::FillIMCombinations(h_IM_2g, 2, photons);


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

    ant::canvas(GetName()+": Basic plots")
            << h_IM_2g
            << h_nPhotons
            << TTree_drawable(t.Tree, "nClusters >> (20,0,20)", "TaggW")
            << endc; // actually draws the canvas
}

AUTO_REGISTER_PHYSICS(Omega_EpEm)
