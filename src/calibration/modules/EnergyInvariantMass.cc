#include "EnergyInvariantMass.h"
#include "analysis/plot/HistogramFactories.h"
#include "analysis/data/Event.h"
#include "analysis/utils/combinatorics.h"

#include "tree/TDataRecord.h"

#include <list>


using namespace std;
using namespace ant;
using namespace ant::calibration;

EnergyInvariantMass::EnergyInvariantMass():
    Calibration::PhysicsModule("EnergyInvariantMass")
{

}

EnergyInvariantMass::ThePhysics::ThePhysics(const string& name):
    Physics(name)
{
    const BinSettings energybins(1000);

    ggIM = HistFac.makeTH1D("ggIM","2 #gamma IM [MeV]","#",energybins,"ggIM");
}


void EnergyInvariantMass::ThePhysics::ProcessEvent(const Event& event)
{
    const auto& cands = event.Reconstructed().Candidates();

    for( auto comb = makeCombination(cands,2); !comb.Done(); ++comb ) {

        if(comb.at(0)->VetoEnergy()==0 && comb.at(1)->VetoEnergy()==0) {
            const Particle a(ParticleTypeDatabase::Photon,comb.at(0));
            const Particle b(ParticleTypeDatabase::Photon,comb.at(1));
            const TLorentzVector gg = a + b;
            ggIM->Fill(gg.M());
        }
    }
}

void EnergyInvariantMass::ThePhysics::Finish()
{
}

void EnergyInvariantMass::ThePhysics::ShowResult()
{
}

unique_ptr<Physics> EnergyInvariantMass::GetPhysicsModule()
{
    return std_ext::make_unique<ThePhysics>(GetName());
}