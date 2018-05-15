#pragma once

#include "calibration/Calibration.h"
#include "base/Detector_t.h"
#include "base/OptionsList.h"
#include "base/Detector_t.h"

#include "tree/TID.h" // for TKeyValue, TID

#include <memory>


namespace ant {

struct ClippedInterpolatorWrapper;

namespace calibration {

class DataManager;

class ClusterCorrection :
        public Calibration::BaseModule,
        public Updateable_traits,
        public ReconstructHook::Clusters
{

public:

    enum class Filter_t {
        MC,
        Data,
        Both
    };

    // ReconstructHook
    virtual void ApplyTo(clusters_t& clusters) override;

    virtual void ApplyTo(TCluster& cluster) =0;

    // Updateable_traits interface
    virtual std::list<Loader_t> GetLoaders() override;

    ClusterCorrection(
            std::shared_ptr<ClusterDetector_t> det,
            const std::string& Name,
            const Filter_t Filter,
            std::shared_ptr<DataManager> calmgr);
    virtual ~ClusterCorrection();

protected:
    const Detector_t::Type_t DetectorType;
    Filter_t filter;

    std::shared_ptr<DataManager> calibrationManager;

    std::unique_ptr<ClippedInterpolatorWrapper> interpolator;

};

/**
 * @brief Cluster energy smearing based on energy and cos(theta)
 */
class ClusterSmearing : public ClusterCorrection {
public:
    using ClusterCorrection::ClusterCorrection;

    void ApplyTo(TCluster& cluster);
};

/**
 * @brief Cluster energy correction based on energy and cluster size
 */
class ClusterECorr : public ClusterCorrection {
public:
    using ClusterCorrection::ClusterCorrection;

    void ApplyTo(TCluster& cluster);
};

class ClusterCorrectionManual : public ClusterCorrection {

public:
    // ReconstructHook
    virtual void ApplyTo(clusters_t& clusters) override;
    virtual void ApplyTo(TCluster& cluster) override =0;

    virtual std::list<Loader_t> GetLoaders() override;

    ClusterCorrectionManual(
            std::shared_ptr<ClusterDetector_t> det,
            const std::string& Name,
            const Filter_t Filter,
            std::shared_ptr<DataManager> calmgr);
    virtual ~ClusterCorrectionManual();
};

/**
 * @brief Manually set additional cluster energy correction factor
 */
class ClusterCorrFactor : public ClusterCorrectionManual {
public:
    ClusterCorrFactor(
            std::shared_ptr<ClusterDetector_t> det,
            const std::string& Name,
            const Filter_t Filter,
            std::shared_ptr<DataManager> calmgr,
            const double corr_factor);

    void ApplyTo(TCluster& cluster);

protected:
    const double factor;
};

/**
 * @brief Manually set cluster energy offset
 */
class ClusterCorrOffset : public ClusterCorrectionManual {
public:

    ClusterCorrOffset(
            std::shared_ptr<ClusterDetector_t> det,
            const std::string& Name,
            const Filter_t Filter,
            std::shared_ptr<DataManager> calmgr,
            const double corr_offset);

    void ApplyTo(TCluster& cluster);

protected:
    const double offset;
};

/**
 * @brief Manually add an additional cluster energy smearing
 */
class ClusterCorrSmearing : public ClusterCorrectionManual {
public:
    ClusterCorrSmearing(
            std::shared_ptr<ClusterDetector_t> det,
            const std::string& Name,
            const Filter_t Filter,
            std::shared_ptr<DataManager> calmgr,
            const double corr_factor);

    void ApplyTo(TCluster& cluster);

protected:
    const double sigma;
};

}}  // namespace ant::calibration
