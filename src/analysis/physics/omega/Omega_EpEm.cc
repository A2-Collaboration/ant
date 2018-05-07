#include "Omega_EpEm.h"
#include "utils/Combinatorics.h"
#include "utils/ParticleTools.h"
#include "base/Logger.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

double Omega_EpEm::effective_radius(const TCandidatePtr& cand) const
{
    return clustertools.EffectiveRadius(*(cand->FindCaloCluster()));
}
double Omega_EpEm::lat_moment(const TCandidatePtr& cand) const
{
    return clustertools.LateralMoment(*(cand->FindCaloCluster()));
}


Omega_EpEm::Omega_EpEm(const string &name, OptionsPtr opts) :
    Physics(name, opts),
    promptrandom(ExpConfig::Setup::Get())
{
    BinSettings bins_nClusters(20);
    BinSettings bins_nParticles(6);
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
    h_nCandCBcharged = HistFac.makeTH1D("charged Candidates/CB",
                                        "Candidates","",
                                        bins_nParticles,
                                        "h_nCandCBcharged"
                                        );
    h_nCandTAPScharged = HistFac.makeTH1D("charged Candidates/TAPS",
                                        "Candidates","",
                                        bins_nParticles,
                                        "h_nCandTAPScharged"
                                        );
    h_nClusters_pr = HistFac.makeTH1D("Number of Clusters - prompt-random",
                                      "nClusters","#",
                                      bins_nClusters,
                                      "h_nClusters_pr"
                                      );
    h_nCombsInIM = HistFac.makeTH1D("Number of Combinations which survive the IM cut",
                                      "nCombsInIM","#",
                                      bins_nParticles,
                                      "h_nCombsInIM"
                                      );
    h_nCandCharged = HistFac.makeTH2D("charged Candidates in TAPS and CB",
                               "charged Candidates in CB",
                               "charged Candidates in TAPS",
                               bins_nParticles,
                               bins_nParticles,
                               "h_nCandCharged"
                               );
    h_nCand = HistFac.makeTH2D("Candidates in TAPS and CB",
                               "Candidates in CB",
                               "Candidates in TAPS",
                               bins_nParticles,
                               bins_nParticles,
                               "h_nCand"
                               );
    h_cbdEE = HistFac.makeTH2D("CB dE-E",
                               "E_{CB} [MeV]",
                               "dE_{PID} [MeV]",
                               BinSettings(1000),
                               BinSettings(100,0,30),
                               "cb_dEE"
                               );
    h_clusteranalysis = HistFac.makeTH2D("Cluster Analysis for candidates in CB",
                               "Lateral Moment",
                               "Effective Radius",
                               BinSettings(100,0,1),
                               BinSettings(100,0,10),
                               "h_clusteranalysis"
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
//        t.Tree->Fill();
    }
    h_nClusters->Fill(event.Reconstructed().Clusters.size()); // how many clusters do we have?

    const auto& recon = event.Reconstructed(); // load reconstructed stuff
    t.IsMC = recon.ID.isSet(TID::Flags_t::MC);
//    const auto& ptree = event.MCTrue().ParticleTree; // load MC true data

    if (t.IsMC) {
        const auto& particletree = event.MCTrue().ParticleTree;
        if (particletree) {
            const auto omegaMC = utils::ParticleTools::FindParticle(ParticleTypeDatabase::Omega, particletree);
            /*const double dE = etapMC->E - sig.etap().E();
            const double theta = std_ext::radian_to_degree(etapMC->Theta());
            h_energy_deviation->Fill(dE);
            h_fsClE_vs_pluto_geant_dE->Fill(dE, etapMC->E - etapMC->M());
            h_theta_vs_vz->Fill(sig.kinfit_ZVertex, theta);
            h_theta_vs_pluto_geant_dE->Fill(dE, theta);
            h_vz_vs_pluto_geant_dE->Fill(dE, sig.kinfit_ZVertex);
            h_delta_vz_vs_pluto_geant_dE->Fill(dE, sig.trueZVertex - sig.kinfit_ZVertex);*/
        } else
            LOG_N_TIMES(1, WARNING) << "(MC debug hists) No particle tree found, only Geant or Pluto file provided, not both";
    }

    if (!triggersimu.HasTriggered()){
//        LOG(INFO) << "Hey, I didn't pass the CBEsum Trigger ";
        hasnttriggered++;
        return; // here we apply the CBEsum Trigger!
    }

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


            }
        }
        else if(ci.Detector & Detector_t::Any_t::CB_Apparatus) {
            cands_cb.emplace_back(i);
            CBSumVetoE += ci.VetoEnergy;
            h_PIDenergy->Fill(ci.VetoEnergy);
            if (ci.VetoEnergy >= 0.3){
                cands_cbCharged.emplace_back(i);

                // now investigate cluster shapes with lateral moment vs effradius
                double effradius = effective_radius(i);
                double latmom = lat_moment(i);
                t.CB_effectiveradius().emplace_back(effradius);
                t.CB_lateralmoment().emplace_back(latmom);
                h_clusteranalysis->Fill(latmom,effradius);
                //LOG(INFO) << "lateral moment: "<< latmom << ", effective radius: " << effradius;
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
    h_nCand->Fill(cands_cb.size(),cands_taps.size());
    h_nCandCBcharged->Fill(cands_cbCharged.size());
    h_nCandTAPScharged->Fill(cands_tapsCharged.size());
    h_nCandCharged->Fill(cands_cbCharged.size(),cands_tapsCharged.size());

    t.nCBneutral = cands_cb.size() - cands_cbCharged.size();
    t.nCBcharged = cands_cbCharged.size();
    t.nTAPSneutral = cands_taps.size() - cands_tapsCharged.size();
    t.nTAPScharged = cands_tapsCharged.size();

    auto m_omega  = ParticleTypeDatabase::Omega.Mass()/1000.0;
    const auto IM_range = interval<double>::CenterWidth(m_omega,100);

    int  nCombsInIM  = 0;
    if (cands_taps.empty()){ // no candidates in TAPS
        if (cands_cbCharged.size() >= 3) { // at least 3 charged candidates in CB

            for( // let's loop over all combinations
                auto comb = analysis::utils::makeCombination(cands_cbCharged,2);
                !comb.done();
                ++comb ) {
                const auto& c1 = comb.at(0);
                const auto& c2 = comb.at(1);
                const TParticle e1(ParticleTypeDatabase::eCharged, c1);
                const TParticle e2(ParticleTypeDatabase::eCharged, c2);
                const auto IM = (e1 + e2).M();
                h_IM->Fill(IM);
                if (IM >= 680 && IM <= 780){
                    nCombsInIM++;
                    t.nCombsInIM++;
                }
                // POOR MAN'S KINFIT: by looking at the mass recoiling against two
                // candidates with the e+e- hypothesis. This should peak at proton mass.
                // const interval<double> MM_cut(850, 1000);
                if(IM_range.Contains(IM)){
                    //
                }
            }
        }
        h_nCombsInIM->Fill(nCombsInIM);
    }
t.fillAndReset(); // do not forget!
}

void Omega_EpEm::ShowResult()
{
    LOG(INFO) << hasnttriggered << "Events didn't trigger";
    ant::canvas(GetName()+": Basic plots")
//            << energy << theta << phi
//            << detectors
//            << h_nClusters
//            << h_nCandidatesEvent
//            << h_nCandCBcharged
//            << h_nCandTAPScharged
//            << h_nCandCB
//            << h_nCandTAPS
//            << drawoption("colz") << h_nCand
//            << drawoption("colz") << h_nCandCharged
//            << h_PIDenergy
//            << h_TAPSVetoEnergy
//            << h_IM
//            << TTree_drawable(t.Tree, "nTAPSneutral >> (8,0,8)")
//            << TTree_drawable(t.Tree, "nTAPScharged >> (8,0,8)")
//            << TTree_drawable(t.Tree, "nCBneutral >> (8,0,8)")
//            << TTree_drawable(t.Tree, "nCBcharged >> (8,0,8)")
            << drawoption("colz") << h_clusteranalysis
//            << drawoption("colz") << TTree_drawable(t.Tree, "nTAPSneutral:nTAPScharged >> (8,0,8,8,0,8)","")
//            << drawoption("colz") << TTree_drawable(t.Tree, "nCBneutral:nCBcharged >> (8,0,8,8,0,8)","")
//            << drawoption("colz") << TTree_drawable(t.Tree, "nTAPSneutral:nCBneutral >> (8,0,8,8,0,8)","")
//            << drawoption("colz") << TTree_drawable(t.Tree, "nTAPScharged:nCBcharged >> (8,0,8,8,0,8)","")
//            << TTree_drawable(t.Tree, "nClusters >> (20,0,20)", "TaggW")
            << endc; // actually draws the canvas
}

AUTO_REGISTER_PHYSICS(Omega_EpEm)
