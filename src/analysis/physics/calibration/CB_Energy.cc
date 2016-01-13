#include "CB_Energy.h"

#include "utils/combinatorics.h"

#include "expconfig/ExpConfig.h"
#include "expconfig/detectors/CB.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::data;
using namespace ant::analysis::physics;

CB_Energy::CB_Energy(const string& name, analysis::PhysOptPtr opts) :
    Physics(name, opts)
{
    auto detector = ExpConfig::Setup::GetDetector(Detector_t::Type_t::CB);

    const BinSettings cb_channels(detector->GetNChannels());
    const BinSettings energybins(1000);

    ggIM = HistFac.makeTH2D("2 neutral IM (CB,CB)", "IM [MeV]", "#",
                            energybins, cb_channels, "ggIM");
    h_cbdisplay = HistFac.make<TH2CB>("h_cbdisplay","Number of entries");
}

void CB_Energy::ProcessEvent(const analysis::data::Event& event)
{
    const auto& cands = event.Reconstructed.Candidates;

    for( auto comb = analysis::utils::makeCombination(cands,2); !comb.Done(); ++comb ) {
        const CandidatePtr& p1 = comb.at(0);
        const CandidatePtr& p2 = comb.at(1);

        if(p1->VetoEnergy==0 && p2->VetoEnergy==0
           && (p1->GetDetector() & Detector_t::Type_t::CB)
           && (p2->GetDetector() & Detector_t::Type_t::CB)) {
            const Particle a(ParticleTypeDatabase::Photon,comb.at(0));
            const Particle b(ParticleTypeDatabase::Photon,comb.at(1));
            const TLorentzVector gg = a + b;

            auto cl1 = p1->FindCaloCluster();
            if(cl1)
                ggIM->Fill(gg.M(),cl1->CentralElement);

            auto cl2 = p2->FindCaloCluster();
            if(cl2)
                ggIM->Fill(gg.M(),cl2->CentralElement);
        }
    }
}

void CB_Energy::ShowResult()
{
    h_cbdisplay->SetElements(*ggIM->ProjectionY());
    canvas(GetName()) << drawoption("colz") << ggIM
                      << h_cbdisplay
                      << endc;
}

AUTO_REGISTER_PHYSICS(CB_Energy)