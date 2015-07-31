#include "TAPSVeto_Energy.h"
#include "analysis/plot/HistogramFactories.h"
#include "analysis/data/Event.h"
#include "analysis/utils/combinatorics.h"

#include "tree/TDataRecord.h"

#include <list>


using namespace std;
using namespace ant;
using namespace ant::calibration;

TAPSVeto_Energy::TAPSVeto_Energy(std::shared_ptr<CalibrationDataManager> calmgr,
                                 Calibration::Converter::ptr_t converter,
                                 double defaultPedestal,
                                 double defaultGain,
                                 double defaultThreshold,
                                 double defaultRelativeGain):
    Energy(Detector_t::Type_t::TAPSVeto,
           calmgr,
           converter,
           defaultPedestal,
           defaultGain,
           defaultThreshold,
           defaultRelativeGain)
{

}

TAPSVeto_Energy::ThePhysics::ThePhysics(const string& name):
    Physics(name)
{
    const BinSettings cb_channels(720);
    const BinSettings energybins(1000);

    ggIM = HistFac.makeTH2D("2 neutral IM (CB,CB)", "IM [MeV]", "#", energybins, cb_channels, "ggIM");
}

void TAPSVeto_Energy::ThePhysics::ProcessEvent(const Event& event)
{
    const auto& cands = event.Reconstructed().Candidates();

    for( auto comb = makeCombination(cands,2); !comb.Done(); ++comb ) {
        const CandidatePtr& p1 = comb.at(0);
        const CandidatePtr& p2 = comb.at(1);

        if(p1->VetoEnergy()==0 && p2->VetoEnergy()==0
           && (p1->Detector() & detector_t::CB)
           && (p2->Detector() & detector_t::CB)) {
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

void TAPSVeto_Energy::ThePhysics::Finish()
{
}

void TAPSVeto_Energy::ThePhysics::ShowResult()
{
    canvas("TAPSVeto_Energy") << drawoption("colz") << ggIM << endc;
}

unique_ptr<Physics> TAPSVeto_Energy::GetPhysicsModule()
{
    return std_ext::make_unique<ThePhysics>(GetName());
}
