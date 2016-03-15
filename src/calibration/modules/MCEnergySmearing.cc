#include "MCEnergySmearing.h"
#include "base/Detector_t.h"
#include "tree/TCluster.h"

#include "TRandom.h"

using namespace std;
using namespace ant;
using namespace ant::calibration;

MCEnergySmearing::MCEnergySmearing(Detector_t::Type_t det, double scale, double smear):
    ReconstructHook::Clusters (), detector(det), scale_factor(scale), smearing(smear)
{}

MCEnergySmearing::~MCEnergySmearing()
{}

void MCEnergySmearing::ApplyTo(ReconstructHook::Base::clusters_t& clusters)
{
    auto myclusters  = clusters.find(detector);

    if(myclusters == clusters.end())
        return;

    for(auto& cluster : myclusters->second) {
        Smear(cluster);
    }
}

void MCEnergySmearing::Smear(TCluster& cluster) const
{
    cluster.Energy = gRandom->Gaus(cluster.Energy*scale_factor, smearing);
}
