#include "Omega_EpEm.h"
#include "utils/Combinatorics.h"
#include "utils/ParticleTools.h"
#include "base/Logger.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;


Omega_EpEm::Omega_EpEm(const string &name, OptionsPtr opts) :
    Physics(name, opts),
    promptrandom(ExpConfig::Setup::Get())
{
    BinSettings bins_nClusters(20);
    BinSettings bins_nParticles(10);
    BinSettings energy_binning(250);
    BinSettings im_binning(200);


    h_nClusters = HistFac.makeTH1D("Number of Clusters", // title
                                   "nClusters","#",      // xlabel, ylabel
                                   bins_nClusters,       // binning
                                   "h_nClusters"         // ROOT object name
                                   );
    h_nCandidatesEvent = HistFac.makeTH1D("Candidates/Event",
                                        "Candidates","",
                                        bins_nClusters,
                                        "nCand"
                                        );
    h_nCandCB = HistFac.makeTH1D("Candidates/CB",
                                        "Candidates","",
                                        bins_nClusters,
                                        "h_nCandCB"
                                        );
    h_nCandTAPS = HistFac.makeTH1D("Candidates/TAPS",
                                        "Candidates","",
                                        bins_nClusters,
                                        "h_nCandTAPS"
                                        );
    h_nClusters_pr = HistFac.makeTH1D("Number of Clusters - prompt-random",
                                      "nClusters","#",
                                      bins_nClusters,
                                      "h_nClusters_pr"
                                      );
    h_cbdEE = HistFac.makeTH2D("CB dE-E",
                               "E_{CB} [MeV]",
                               "dE_{PID} [MeV]",
                               BinSettings(1000),
                               BinSettings(100,0,30),
                               "cb_dEE"
                               );
    h_PIDenergy = HistFac.makeTH1D("Energy in PID","E [MeV]","",BinSettings(100,0,10),"h_PIDenergy");
    h_TAPSVetoEnergy = HistFac.makeTH1D("Energy in TAPS Veto","E [MeV]","",BinSettings(100,0,10),"h_TAPSVetoEnergy");
    h_IM = HistFac.makeTH1D("IM of two charged particles","E [MeV]","",BinSettings(1000,0,1000),"h_IM");
    energy = HistFac.makeTH1D("Energy","E [MeV]","",BinSettings(1000),"energy");
    theta  = HistFac.makeTH1D("Theta","#theta [#circ]","",BinSettings(360,0,180),"theta");
    phi    = HistFac.makeTH1D("Phi","#phi [#circ]","",BinSettings(720,-180,180),"phi");
    detectors = HistFac.makeTH1D("Detectors","","", BinSettings(1),"detectors");

    // prompt random window
    //promptrandom.AddPromptRange({ -7,   7}); // in ns
    //promptrandom.AddRandomRange({-50, -10});
    //promptrandom.AddRandomRange({ 10,  50});

    // some tree
    t.CreateBranches(HistFac.makeTTree("t"));

    LOG(INFO) << "Initialized " << GetName() << "";
}

void Omega_EpEm::ProcessEvent(const TEvent& event, manager_t&)
{
    if(!triggersimu.ProcessEvent(event)) {
        return;
    }
    for(auto& taggerhit : event.Reconstructed().TaggerHits) {
        promptrandom.SetTaggerTime(triggersimu.GetCorrectedTaggerTime(taggerhit));
        if(promptrandom.State() == PromptRandom::Case::Outside) {
            continue;
        }
        h_nClusters_pr->Fill(event.Reconstructed().Clusters.size(), promptrandom.FillWeight());

        t.TaggW = promptrandom.FillWeight();
        t.nClusters = event.Reconstructed().Clusters.size();
        t.Tree->Fill();
    }
    h_nClusters->Fill(event.Reconstructed().Clusters.size()); // how many clusters do we have?


    const auto& recon = event.Reconstructed();
    t.IsMC = recon.ID.isSet(TID::Flags_t::MC);

    // get candidate list and make lists for candidates in CB and TAPS apparatus
    const auto& cands = event.Reconstructed().Candidates;

    TCandidatePtrList cands_taps;
    TCandidatePtrList cands_tapsCharged;
    TCandidatePtrList cands_cb;
    TCandidatePtrList cands_cbCharged;

    h_nCandidatesEvent->Fill(cands.size()); // how many candidates do we have?
    CBSumVetoE = 0;

    auto i = cands.begin();
    while(i!=cands.end()) {
        const TCandidate& ci = *i; // every candidate is now ci, the candidate pointer is i

        if(ci.Detector & Detector_t::Any_t::TAPS_Apparatus) {
            cands_taps.emplace_back(i); // using the pointer because this is a pointerlist
            h_TAPSVetoEnergy->Fill(ci.VetoEnergy);
            if (ci.VetoEnergy >= 0.5){
                cands_tapsCharged.emplace_back(i);
                // kann ich jetzt hier ci in den tree packen? hmmmmm
//                t.p_tapsCharged().emplace_back(*i);

            }
        }
        else if(ci.Detector & Detector_t::Any_t::CB_Apparatus) {
            cands_cb.emplace_back(i);
            CBSumVetoE += ci.VetoEnergy;
            h_PIDenergy->Fill(ci.VetoEnergy);
            if (ci.VetoEnergy >= 0.3){
                cands_cbCharged.emplace_back(i);

            }
        }
        energy->Fill(ci.CaloEnergy); // deponierte Energie in den Kalorimetern
        theta->Fill(std_ext::radian_to_degree(ci.Theta));
        phi->Fill(std_ext::radian_to_degree(ci.Phi));
        detectors->Fill(string(ci.Detector).c_str(),1);

        ++i;
    }

    h_nCandCB->Fill(cands_cb.size());
    h_nCandTAPS->Fill(cands_taps.size());

//    Int_t nCBneutral, nTAPSneutral;
//    nCBneutral = cands_cbCharged.size() - cands_cb.size();
//    nTAPSneutral = cands_tapsCharged.size() - cands_taps.size();

    if (cands_taps.empty()){
        if (cands_cbCharged.size() >= 3) {

            for(
                auto comb = analysis::utils::makeCombination(cands_cbCharged,2);
                !comb.done();
                ++comb ) {
                const auto& c1 = comb.at(0);
                const auto& c2 = comb.at(1);
                const TParticle e1(ParticleTypeDatabase::eCharged, c1);
                const TParticle e2(ParticleTypeDatabase::eCharged, c2);
                h_IM->Fill((e1 + e2).M());
            }
        }
    }

}

void Omega_EpEm::ShowResult()
{

    ant::canvas(GetName()+": Basic plots")
            << energy << theta << phi
            << detectors
//            << h_nClusters
            << h_nCandidatesEvent
            << h_nCandCB
            << h_nCandTAPS
            << h_PIDenergy
            << h_TAPSVetoEnergy
            << h_IM
            << TTree_drawable(t.Tree, "nClusters >> (20,0,20)", "TaggW")
            << endc; // actually draws the canvas
}

AUTO_REGISTER_PHYSICS(Omega_EpEm)
