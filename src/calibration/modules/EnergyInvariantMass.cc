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
    Calibration::Module("EnergyInvariantMass")
{
    const BinSettings energybins(1000);

    ggIM = HistFac.makeTH1D("ggIM","2 #gamma IM [MeV]","#",energybins,"ggIM");
}





void EnergyInvariantMass::ProcessEvent(const Event &event)
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

void EnergyInvariantMass::Finish()
{

}

void EnergyInvariantMass::ShowResult()
{

}

void EnergyInvariantMass::ApplyTo(const readhits_t& hits, extrahits_t&)
{

}

std::list<TID> EnergyInvariantMass::GetChangePoints() const
{
    return {};
}

void EnergyInvariantMass::Update(const TID &id)
{

}

AUTO_REGISTER_PHYSICS(EnergyInvariantMass, "EnergyInvariantMass")
