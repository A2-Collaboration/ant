#pragma once

#include "calibration/Calibration.h"
#include "base/Detector_t.h"
#include "base/OptionsList.h"
#include "base/Detector_t.h"

#include "tree/TID.h" // for TKeyValue, TID

#include <memory>


namespace ant {

namespace calibration {

class DataManager;

class ClusterCorrection :
        public Calibration::BaseModule,
        public Updateable_traits,
        public ReconstructHook::Clusters
{

public:
    // ReconstructHook
    virtual void ApplyTo(clusters_t& clusters) override;

    virtual void ApplyTo(TCluster& cluster) =0;

    // Updateable_traits interface
    virtual std::list<Loader_t> GetLoaders() override;

    ClusterCorrection(
            std::shared_ptr<ClusterDetector_t> det,
            const std::string& Name,
           std::shared_ptr<DataManager> calmgr);
    virtual ~ClusterCorrection();

protected:
    const Detector_t::Type_t DetectorType;

    std::shared_ptr<DataManager> calibrationManager;

    struct Interpolator;
    std::unique_ptr<Interpolator> interpolator;

};

class ClusterSmearing : public ClusterCorrection {
public:
    using ClusterCorrection::ClusterCorrection;

    void ApplyTo(TCluster& cluster);
};

class ClusterScaling : public ClusterCorrection {
public:
    using ClusterCorrection::ClusterCorrection;

    void ApplyTo(TCluster& cluster);
};

}}  // namespace ant::calibration
