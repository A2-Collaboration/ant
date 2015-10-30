#pragma once

#include "base/Detector_t.h"

#include "tree/MemoryPool.h"
#include "base/mapped_vectors.h"

#include <memory>
#include <map>
#include <list>

namespace ant {

class TID;
class THeaderInfo;
class TDetectorRead;
class TDetectorReadHit;
class TEvent;
class TCluster;

namespace reconstruct {
class AdaptorTClusterHit;
}

struct Reconstruct_traits {
    virtual ~Reconstruct_traits() = default;
    /**
     * @brief Initialize called before first DoReconstruct() call
     * @param headerInfo
     */
    virtual void Initialize(const THeaderInfo& headerInfo) = 0;
    /**
     * @brief DoReconstruct shall convert the given TDetectorRead to TEvent
     * @param detectorRead
     * @return the reconstructed TEvent
     */
    virtual MemoryPool<TEvent>::Item DoReconstruct(TDetectorRead& detectorRead) = 0;
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
     * @brief The Base class just defines some useful types
     */
    class Base {
    public:
        using readhits_t = std_ext::mapped_vectors< Detector_t::Type_t, TDetectorReadHit* >;
        using extrahits_t = std::list< TDetectorReadHit >;
        using clusterhits_t = std::map< Detector_t::Type_t, std::list< reconstruct::AdaptorTClusterHit > >;
        using clusters_t = std::map< Detector_t::Type_t, std::list< TCluster > >;
        virtual ~Base() = default;
    };

    /**
     * @brief The DetectorReadHits class instances are applied before hit matching
     */
    class DetectorReadHits : public Base {
    public:
        virtual void ApplyTo(const readhits_t& hits, extrahits_t& extrahits) = 0;
    };

    /**
     * @brief The ClusterHits class instances are applied before clustering and after hit matching
     */
    class ClusterHits : public Base {
    public:
        virtual void ApplyTo(clusterhits_t& clusterhits) = 0;
    };

    /**
     * @brief The Clusters class instances are applied before candidate matching and after clustering
     */
    class Clusters : public Base {
    public:
        virtual void ApplyTo(clusters_t& clusters) = 0;
    };
};


/**
 * @brief The Updateable_traits class used by the UpdateableManager
 */
class Updateable_traits {
public:
    /**
     * @brief GetChangePoints builds (possibly unsorted) list of points
     * when Updateable_traits::Update() should be called
     * @return vectorized list of change points, index is given to Update()
     */
    virtual std::vector<std::list<TID>> GetChangePoints() const = 0;
    /**
     * @brief Update shall load new parameters for the upcoming event ids
     * @param index of element in vector returned by GetChangePoints() or UpdateOnFirstEvent()
     * @param id indicates moment of change
     */
    virtual void Update(std::size_t index, const TID& id) = 0;
    /**
     * @brief UpdatedTIDFlags called when processed event has some different flags in TID
     * @param id the ID with some different Flags field
     */
    virtual void UpdatedTIDFlags(const TID&) {}
};

} // namespace ant
