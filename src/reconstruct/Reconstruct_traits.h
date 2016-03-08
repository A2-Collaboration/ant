#pragma once

#include "base/Detector_t.h"
#include "base/mapped_vectors.h"

#include <memory>
#include <map>
#include <list>
#include <functional>

namespace ant {

struct TID;
struct TEventData;
struct TDetectorReadHit;

struct TClusterHit;
using TClusterHitList = std::vector<TClusterHit>;

struct TCluster;
using TClusterPtr = std::shared_ptr<TCluster>;
using TClusterList = std::vector<TClusterPtr>;

struct Reconstruct_traits {
    /**
     * @brief DoReconstruct shall convert the given TDetectorRead to TEvent
     * @param detectorRead
     * @return the reconstructed TEvent
     */
    virtual void DoReconstruct(TEventData& reconstructed) = 0;

    virtual ~Reconstruct_traits() = default;
};

/**
 * @brief The ReconstructHook struct
 *
 * A reconstruct hook is a calibration most of the time, which modifies
 * the information in the detector read hits from raw data to meaningful
 * physical quantities. However, the concept is more general and also regards
 * for example the TAPS Shower Correction as a reconstruct hook for clusters,
 * or some of the calibration/converters hook into the reconstruct in order to
 * find their reference timings.
 *
 */
struct ReconstructHook {
    /**
     * @brief The Base struct just defines some useful types
     */
    struct Base {
        using readhits_t = std_ext::mapped_vectors< Detector_t::Type_t, TDetectorReadHit* >;
        using clusterhits_t = std::map< Detector_t::Type_t, TClusterHitList >;
        using clusters_t = std::map< Detector_t::Type_t, TClusterList >;
        virtual ~Base() = default;
    };

    /**
     * @brief The DetectorReadHits struct instances are applied before hit matching
     */
    struct DetectorReadHits : virtual Base {
        virtual void ApplyTo(const readhits_t& hits) = 0;
    };

    /**
     * @brief The ClusterHits struct instances are applied before clustering and after hit matching
     */
    struct ClusterHits : virtual Base {
        virtual void ApplyTo(clusterhits_t& clusterhits) = 0;
    };

    /**
     * @brief The Clusters struct instances are applied before candidate matching and after clustering
     */
    struct Clusters : virtual Base {
        virtual void ApplyTo(clusters_t& clusters) = 0;
    };

    /**
     * @brief The EventData struct instances are applied after candidate building
     */
    struct EventData : virtual Base {
        virtual void ApplyTo(TEventData& reconstructed) = 0;
    };

};


/**
 * @brief The Updateable_traits class used by the UpdateableManager
 */
struct Updateable_traits {

    using Loader_t = std::function<void(const TID& currPoint, TID& nextChangePoint)>;

    virtual std::list<Loader_t> GetLoaders() = 0;

    /**
     * @brief UpdatedTIDFlags called when processed event has some different flags in TID
     * @param id the ID with some different Flags field
     */
    virtual void UpdatedTIDFlags(const TID&) {}

    virtual ~Updateable_traits() = default;
};

} // namespace ant
