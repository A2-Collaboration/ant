#include "TAPS_Energy.h"

#include "utils/combinatorics.h"

#include "base/std_ext/vector.h"

#include "expconfig/ExpConfig.h"
#include "expconfig/detectors/TAPS.h"

#include "TTree.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::data;
using namespace ant::analysis::physics;



TAPS_Energy::TAPS_Energy(const string& name, analysis::PhysOptPtr opts) :
    Physics(name, opts)
{
    taps_detector = ExpConfig::Setup::GetDetector<expconfig::detector::TAPS>();

    const BinSettings taps_channels(taps_detector->GetNChannels());
    const BinSettings energybins(1000);
    const BinSettings timebins(1000,-100,100);

    ggIM = HistFac.makeTH2D("2 neutral IM (TAPS,CB)", "IM [MeV]", "#",
                            energybins, taps_channels, "ggIM");


    timing_cuts = HistFac.makeTH2D("Check timing cuts", "Time [ns]", "#",
                                   timebins, taps_channels, "timing_cuts");

    h_pedestals = HistFac.makeTH2D(
                      "TAPS Pedestals",
                      "Raw ADC value",
                      "#",
                      BinSettings(300),
                      taps_channels,
                      "Pedestals");

    if(!opts->Get<bool>("TAPS_Energy_Advanced"))
        return;

    ggIM_mult = HistFac.makeTH3D("IM (TAPS,CB) mult",
                                 "IM [MeV]", "Channel", "Multiplicity",
                                 energybins, taps_channels, BinSettings(10),
                                 "ggIM_mult");

    cands_tree = HistFac.makeTTree("cands_tree");
    cands_CB.Setup("CB", cands_tree);
    cands_TAPS.Setup("TAPS", cands_tree);
}

void TAPS_Energy::ProcessEvent(const Event& event)
{
    const auto& cands = event.Reconstructed.Candidates;

    // pedestals
    for(const Cluster& cluster : event.Reconstructed.AllClusters) {
        if(!(cluster.Detector & Detector_t::Type_t::TAPS))
            continue;
        for(const Cluster::Hit& clusterhit : cluster.Hits) {
            /// \todo check for timing hit?
            /// \todo check for trigger pattern?
            for(const Cluster::Hit::Datum& datum : clusterhit.Data) {
                if(datum.Type != Channel_t::Type_t::Pedestal)
                    continue;
                h_pedestals->Fill(datum.Value, clusterhit.Channel);
            }
        }
    }

    // invariant mass of two photons
    for( auto comb = analysis::utils::makeCombination(cands,2); !comb.Done(); ++comb ) {
        const CandidatePtr& cand1 = comb.at(0);
        const CandidatePtr& cand2 = comb.at(1);

        if(cand1->VetoEnergy<0.5 && cand2->VetoEnergy<0.5) {

            //require exactly 1 CB and 1 TAPS
            const auto CBTAPS = Detector_t::Type_t::CB | Detector_t::Type_t::TAPS;
            const auto dets = (cand1->Detector & CBTAPS) ^ (cand2->Detector & CBTAPS);

            if(dets & CBTAPS) {
                const Particle a(ParticleTypeDatabase::Photon,comb.at(0));
                const Particle b(ParticleTypeDatabase::Photon,comb.at(1));
                const TLorentzVector& gg = a + b;


                // Find the one that was in TAPS
                auto cand_taps = cand1->Detector & Detector_t::Type_t::TAPS ? cand1 : cand2;
                auto cl_taps = cand_taps->FindCaloCluster();
                if(cl_taps) {
                    const unsigned ch = cl_taps->CentralElement;
                    const unsigned ring = taps_detector->GetRing(ch);

                    // fill in IM only if ring>4 or if timecut is passed
                    double weight = -1.0;
                    if(ring > 4 || fabs(cand_taps->Time) < 5) {
                        weight = 1.0;
                        ggIM->Fill(gg.M(),ch);
                    }
                    timing_cuts->Fill(cand_taps->Time, ch, weight);
                }
            }
        }
    }

    if(ggIM_mult != nullptr) {
        // collect neutral candidates
        CandidateList cb_candidates;
        CandidateList taps_candidates;

        for(const CandidatePtr& cand : event.Reconstructed.Candidates)
        {
            if((cand->Detector & Detector_t::Type_t::TAPS)
               && cand->VetoEnergy < 1.0)
            {
                auto taps_cluster = cand->FindCaloCluster();
                const unsigned ch = taps_cluster->CentralElement;
                const unsigned ring = taps_detector->GetRing(ch);
                if(ring > 4 || fabs(cand->Time) < 2.5)
                    taps_candidates.emplace_back(move(cand));
            }
            else if((cand->Detector & Detector_t::Type_t::CB)
                    && cand->VetoEnergy < 0.25)
            {
                cb_candidates.emplace_back(move(cand));
            }
        }

        const auto nNeutrals = taps_candidates.size() + cb_candidates.size();

        for(const CandidatePtr& taps_cand : taps_candidates) {
            auto taps_cluster = taps_cand->FindCaloCluster();
            const unsigned ch = taps_cluster->CentralElement;
            for(const CandidatePtr& cb_cand : cb_candidates) {

                const TLorentzVector& gg = Particle(ParticleTypeDatabase::Photon, taps_cand)
                                           + Particle(ParticleTypeDatabase::Photon, cb_cand);
                ggIM_mult->Fill(gg.M(), ch, nNeutrals);
            }
        }
    }

    if(cands_tree != nullptr) {
        // fill some tree
        cands_TAPS.Clear();
        cands_CB.Clear();
        bool interesting_event = false;
        const vector<unsigned> interesting_channels = {31, 0, 1, 2, 15};
        for(const CandidatePtr& cand : event.Reconstructed.Candidates) {
            // find interesting TAPS clusters
            if(cand->Detector & Detector_t::Type_t::TAPS) {
                auto taps_cluster = cand->FindCaloCluster();
                if(std_ext::contains(interesting_channels, taps_cluster->CentralElement)) {
                    cands_TAPS.Fill(*cand);
                    interesting_event = true;
                }
            }
        }

        if(interesting_event) {
            bool CB_seen = false;
            for(const CandidatePtr& cand : event.Reconstructed.Candidates) {
                if(cand->Detector & Detector_t::Type_t::CB)
                {
                    cands_CB.Fill(*cand);
                    CB_seen = true;
                }
            }
            if(CB_seen)
                cands_tree->Fill();
        }
    }
}

void TAPS_Energy::ShowResult()
{
    canvas(GetName()) << drawoption("colz") << ggIM
                      << drawoption("colz") << timing_cuts
                      << drawoption("colz") << h_pedestals
                      << endc;
    if(ggIM_mult) {
        canvas c_proj(GetName()+" Multiplicities");
        for(unsigned n=2;n<=10;n++) {
            ggIM_mult->GetZaxis()->SetRange(n,n+1);
            stringstream ss_name;
            ss_name << "Mult_" << n << "_yx";
            TH2* proj = dynamic_cast<TH2*>(ggIM_mult->Project3D(ss_name.str().c_str()));
            c_proj << drawoption("colz") << proj;
        }
        c_proj << endc;
    }
}



void TAPS_Energy::tree_data_t::Setup(const string& prefix, TTree* tree)
{
    tree->Branch((prefix + "_Ek").c_str(),      &Ek);
    tree->Branch((prefix + "_Theta").c_str(),   &Theta);
    tree->Branch((prefix + "_Phi").c_str(),     &Phi);
    tree->Branch((prefix + "_VetoE").c_str(),   &VetoE);
    tree->Branch((prefix + "_Time").c_str(),    &Time);
    tree->Branch((prefix + "_Channel").c_str(), &Channel);

}

void TAPS_Energy::tree_data_t::Clear()
{
    Ek.resize(0);
    Theta.resize(0);
    Phi.resize(0);
    VetoE.resize(0);
    Time.resize(0);
    Channel.resize(0);
}

void TAPS_Energy::tree_data_t::Fill(const Candidate& cand)
{
    Ek.push_back(cand.ClusterEnergy);
    Theta.push_back(cand.Theta);
    Phi.push_back(cand.Phi);
    VetoE.push_back(cand.VetoEnergy);
    Time.push_back(cand.Time);
    auto cl = cand.FindCaloCluster();
    Channel.push_back(cl->CentralElement);
}

AUTO_REGISTER_PHYSICS(TAPS_Energy)
