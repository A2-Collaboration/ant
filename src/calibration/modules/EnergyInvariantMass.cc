#include "EnergyInvariantMass.h"
#include "analysis/plot/HistogramFactories.h"
#include "analysis/data/Event.h"
#include "analysis/utils/combinatorics.h"

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
    const auto& photons = event.Reconstructed().Particles().Get(ParticleTypeDatabase::Photon);

    for( auto comb = makeCombination(photons,2); !comb.Done(); ++comb ) {

        TLorentzVector gg = *comb.at(0) + *comb.at(1);
        ggIM->Fill(gg.M());
    }
}

void EnergyInvariantMass::Finish()
{

}

void EnergyInvariantMass::ShowResult()
{

}

void EnergyInvariantMass::ApplyTo(const std::map< Detector_t::Type_t, std::list< TDetectorReadHit* > >& hits)
{

}

void EnergyInvariantMass::BuildRanges(list<TID>& ranges)
{

}

void EnergyInvariantMass::Update(const TID &id)
{

}