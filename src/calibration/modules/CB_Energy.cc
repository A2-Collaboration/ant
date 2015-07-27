#include "CB_Energy.h"
#include "analysis/plot/HistogramFactories.h"
#include "analysis/data/Event.h"
#include "analysis/utils/combinatorics.h"

#include "tree/TDataRecord.h"

#include <list>


using namespace std;
using namespace ant;
using namespace ant::calibration;

CB_Energy::CB_Energy(Calibration::Converter::ptr_t converter,
                               const double defaultPedestal,
                               const double defaultGain,
                               const double defaultThreshold):
    Integral(Detector_t::Type_t::CB,
             converter,
             defaultPedestal,
             defaultGain,
             defaultThreshold)
{

}

CB_Energy::ThePhysics::ThePhysics(const string& name):
    Physics(name)
{
    const BinSettings energybins(1000);

    ggIM = HistFac.makeTH1D("ggIM","2 #gamma IM [MeV]","#",energybins,"ggIM");
}


void CB_Energy::ThePhysics::ProcessEvent(const Event& event)
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
            ggIM->Fill(gg.M());
        }
    }
}

void CB_Energy::ThePhysics::Finish()
{
}

void CB_Energy::ThePhysics::ShowResult()
{
}

unique_ptr<Physics> CB_Energy::GetPhysicsModule()
{
    return std_ext::make_unique<ThePhysics>(GetName());
}
