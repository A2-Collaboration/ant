#include "CB_Energy.h"

#include "utils/Combinatorics.h"

#include "expconfig/ExpConfig.h"
#include "expconfig/detectors/CB.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::physics;

void CB_Energy::FillggIM(const TCluster &cl1, const TCluster &cl2, const double imass)
{
    if(!RequireClean || (!cl2.HasFlag(TCluster::Flags_t::TouchesHoleCentral)))
        ggIM->Fill(imass, cl1.CentralElement);
}

CB_Energy::CB_Energy(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    RequireClean(opts->Get<bool>("RequireClean", true)),
    RequireVetoEZero(opts->Get<bool>("RequireVetoEZero", true)),
    MinOpeningAngle(opts->Get<double>("MinOpeningAngle", 0))

{
    auto detector = ExpConfig::Setup::GetDetector(Detector_t::Type_t::CB);

    const BinSettings cb_channels(detector->GetNChannels());
    const BinSettings energybins(1000);

    ggIM = HistFac.makeTH2D("2 neutral IM (CB,CB)", "IM [MeV]", "#",
                            energybins, cb_channels, "ggIM");
    h_cbdisplay = HistFac.make<TH2CB>("h_cbdisplay","Number of entries");
}

void CB_Energy::ProcessEvent(const TEvent& event, manager_t&)
{
    const auto& cands = event.Reconstructed().Candidates;

    for( auto comb = analysis::utils::makeCombination(cands.get_ptr_list(),2); !comb.done(); ++comb ) {
        const TCandidatePtr& p1 = comb.at(0);
        const TCandidatePtr& p2 = comb.at(1);

        if(   (!RequireVetoEZero || p1->VetoEnergy==0)
           && (!RequireVetoEZero || p2->VetoEnergy==0)
           && (p1->Detector & Detector_t::Type_t::CB)
           && (p2->Detector & Detector_t::Type_t::CB)
           && (vec3(*p1).Angle(*p2) > MinOpeningAngle)
           )
        {
            const TParticle g1(ParticleTypeDatabase::Photon,p1);
            const TParticle g2(ParticleTypeDatabase::Photon,p2);
            const auto ggmass = (g1 + g2).M();

            const auto cl1 = p1->FindCaloCluster();
            const auto cl2 = p2->FindCaloCluster();

            if(cl1 && cl2) {
                FillggIM(*cl1, *cl2, ggmass);
                FillggIM(*cl2, *cl1, ggmass);
            }
        }
    }
}

void CB_Energy::ShowResult()
{
    auto proj = dynamic_cast<TH1D*>(ggIM->ProjectionX());
    proj->GetXaxis()->SetRangeUser(0, 300);
    h_cbdisplay->SetElements(*ggIM->ProjectionY());
    canvas(GetName()) << drawoption("colz") << ggIM
                      << h_cbdisplay
                      << proj
                      << endc;
}

AUTO_REGISTER_PHYSICS(CB_Energy)
