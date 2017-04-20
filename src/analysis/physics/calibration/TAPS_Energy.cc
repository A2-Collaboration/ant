#include "TAPS_Energy.h"

#include "utils/Combinatorics.h"

#include "base/std_ext/container.h"

#include "expconfig/ExpConfig.h"
#include "expconfig/detectors/TAPS.h"


#include "TTree.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::physics;



TAPS_Energy::TAPS_Energy(const string& name, OptionsPtr opts) :
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

    h_tapsdisplay = HistFac.make<TH2TAPS>("h_tapsdisplay","Number of entries");


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

void TAPS_Energy::ProcessEvent(const TEvent& event, manager_t&)
{
    const auto& cands = event.Reconstructed().Candidates;

    // pedestals
    for(const TDetectorReadHit& readhit : event.Reconstructed().DetectorReadHits) {
        if(readhit.DetectorType != Detector_t::Type_t::TAPS)
            continue;
        if(readhit.ChannelType != Channel_t::Type_t::Integral)
            continue;
        /// \todo check for timing hit?
        /// \todo check for trigger pattern?
        for(const auto& value : readhit.Values)
            h_pedestals->Fill(value.Uncalibrated, readhit.Channel);
    }

    // invariant mass of two photons
    for( auto comb = analysis::utils::makeCombination(cands.get_ptr_list(),2); !comb.done(); ++comb ) {
        const TCandidatePtr& p1 = comb.at(0);
        const TCandidatePtr& p2 = comb.at(1);

        // Use PID for CB, but not Vetos for TAPS
        const auto checkVeto = [] (const TCandidate& c) {
            return (c.Detector & Detector_t::Type_t::TAPS) || (c.VetoEnergy < 0.5);
        };

        if(checkVeto(*p1) && checkVeto(*p2)) {

            {

                const TParticle g1(ParticleTypeDatabase::Photon,comb.at(0));
                const TParticle g2(ParticleTypeDatabase::Photon,comb.at(1));
                const auto& gg = (g1 + g2).M();

                const auto Fill = [&gg,this] (const TCandidate& c1, const TCandidate& c2) {

                    if(!(c1.Detector & Detector_t::Type_t::TAPS))
                        return;

                    auto cl1 = c1.FindCaloCluster();
                    auto cl2 = c2.FindCaloCluster();

                    if(cl1 && cl2 && !cl2->HasFlag(TCluster::Flags_t::TouchesHoleCentral)) {
                        const unsigned ch = cl1->CentralElement;
                        const unsigned ring = taps_detector->GetRing(ch);

                        // fill in IM only if ring>4 or if timecut is passed
                        double weight = -1.0;
                        if(ring > 4 || fabs(cl1->Time) < 5) {
                            weight = 1.0;
                            ggIM->Fill(gg,ch);
                        }
                        timing_cuts->Fill(cl1->Time, ch, weight);
                    }
                };

                Fill(*p1,*p2);
                Fill(*p2,*p1);

            }
        }
    }

    if(ggIM_mult != nullptr) {
        // collect neutral candidates
        TCandidatePtrList cb_candidates;
        TCandidatePtrList taps_candidates;

        for(const auto& cand : event.Reconstructed().Candidates.get_iter())
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

        for(const TCandidatePtr& taps_cand : taps_candidates) {
            auto taps_cluster = taps_cand->FindCaloCluster();
            const unsigned ch = taps_cluster->CentralElement;
            for(const TCandidatePtr& cb_cand : cb_candidates) {

                const auto& gg = TParticle(ParticleTypeDatabase::Photon, taps_cand)
                                 + TParticle(ParticleTypeDatabase::Photon, cb_cand);
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
        for(const TCandidate& cand : event.Reconstructed().Candidates) {
            // find interesting TAPS clusters
            if(cand.Detector & Detector_t::Type_t::TAPS) {
                auto taps_cluster = cand.FindCaloCluster();
                if(std_ext::contains(interesting_channels, taps_cluster->CentralElement)) {
                    cands_TAPS.Fill(cand);
                    interesting_event = true;
                }
            }
        }

        if(interesting_event) {
            bool CB_seen = false;
            for(const TCandidate& cand : event.Reconstructed().Candidates) {
                if(cand.Detector & Detector_t::Type_t::CB)
                {
                    cands_CB.Fill(cand);
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
    h_tapsdisplay->SetElements(*ggIM->ProjectionY());

    canvas(GetName()) << drawoption("colz") << ggIM
                      << drawoption("colz") << timing_cuts
                      << drawoption("colz") << h_pedestals
                      << h_tapsdisplay
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

void TAPS_Energy::tree_data_t::Fill(const TCandidate& cand)
{
    Ek.push_back(cand.CaloEnergy);
    Theta.push_back(cand.Theta);
    Phi.push_back(cand.Phi);
    VetoE.push_back(cand.VetoEnergy);
    Time.push_back(cand.Time);
    auto cl = cand.FindCaloCluster();
    Channel.push_back(cl->CentralElement);
}

AUTO_REGISTER_PHYSICS(TAPS_Energy)
