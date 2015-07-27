#pragma once

#include "calibration/Calibration.h"
#include "tree/TCluster.h"
#include <cmath>

namespace ant {
namespace calibration {

class TAPS_ShowerCorrection :
        public Calibration::BaseModule,
        public ReconstructHook::Clusters
{
public:
    TAPS_ShowerCorrection() : Calibration::BaseModule("TAPS_ShowerCorrection") {}

    virtual void ApplyTo(clusters_t& sorted_clusters) override {
        // search for TAPS clusters
        const auto it_sorted_clusters = sorted_clusters.find(Detector_t::Type_t::TAPS);
        if(it_sorted_clusters == sorted_clusters.end())
            return;

        std::list< TCluster >& clusters = it_sorted_clusters->second;

        for(TCluster& cluster : clusters) {
            const double sh_dep = 2.05 * (std::log(cluster.Energy / 12.7) + 1.2);
            if (sh_dep > 0)
            {
                TVector3& pos = cluster.Position;
                const double sh_corr = pos.Mag() / sh_dep + 1.;
                pos.SetX(pos.X() - pos.X() / sh_corr);
                pos.SetY(pos.Y() - pos.Y() / sh_corr);
            }
        }
    }
};


}}