#include "Uncertainties.h"

#include "expconfig/ExpConfig.h"
#include "base/Logger.h"

using namespace std;
using namespace ant::analysis::utils;

UncertaintyModel::UncertaintyModel()
{
    try {
        tagger = ExpConfig::Setup::GetDetector<TaggerDetector_t>();
    }
    catch(ExpConfig::Exception&) {
        LOG(WARNING) << "No tagger found in setup, using default uncertainties (probably NOT what you want!)";
    }
}

UncertaintyModel::~UncertaintyModel()
{}

double UncertaintyModel::GetBeamEnergySigma(double photon_energy) const
{
    if(tagger) {
        unsigned channel;
        if(tagger->TryGetChannelFromPhoton(photon_energy, channel)) {
            // sigma of uniform distribution is width/sqrt(12)
            return tagger->GetPhotonEnergyWidth(channel)/sqrt(12.0);
        }
    }

    // default uncertainty 1MeV
    return 1;
}





