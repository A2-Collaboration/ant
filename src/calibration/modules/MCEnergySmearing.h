#pragma once

#include "base/Detector_t.h"
#include "reconstruct/Reconstruct_traits.h"

namespace ant {
namespace calibration {

/**
 * @brief Basic Energy smearing for MC data. Works on clusters and scales the energy and then smears it with a gaussian.
 */
class MCEnergySmearing : public ReconstructHook::Clusters {
protected:
    const ant::Detector_t::Type_t detector;
    const double scale_factor;
    const double smearing;

public:
    /**
     * @brief MCEnergySmearing
     * @param det The detector type to work on
     * @param scale scale factor for the energy
     * @param smear the sigma of the gaussian to smear with
     */
    MCEnergySmearing(ant::Detector_t::Type_t det, double scale, double smear);
    virtual ~MCEnergySmearing();

    virtual void ApplyTo(clusters_t& clusters) override;

    /**
     * @brief Smearinf function. Modifies the cluster energy.
     * @param cluster
     *
     * This function is a conveninet place to override and implement your own smearing function.
     */
    virtual void Smear(TCluster& cluster) const;
};

}
}
