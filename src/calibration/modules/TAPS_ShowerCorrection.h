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

        auto& clusters = it_sorted_clusters->second;

        for(TCluster& cluster : clusters) {
            const auto sh_dep = 2.05 * (std::log(cluster.Energy / 12.7) + 1.2);
            if (sh_dep > 0)
            {
                auto& pos = cluster.Position;
                const auto sh_corr = pos.R() / sh_dep + 1.;
                pos.x -= pos.x / sh_corr;
                pos.y -= pos.y / sh_corr;
            }
        }
    }
};


}}