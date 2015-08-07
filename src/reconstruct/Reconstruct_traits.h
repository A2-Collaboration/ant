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


struct Reconstruct_traits {
    virtual ~Reconstruct_traits() = default;
    virtual void Initialize(const THeaderInfo& headerInfo) = 0;
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
 * find their reference timings or scalers.
 *
 */
struct ReconstructHook {

    class Base {
    public:
//        using readhits_t = std::map< Detector_t::Type_t, std::list< TDetectorReadHit* > >;
        using readhits_t = std_ext::mapped_vectors< Detector_t::Type_t, TDetectorReadHit* >;
        using extrahits_t = std::list< TDetectorReadHit >;
        using clusters_t = std::map< Detector_t::Type_t, std::list< TCluster > >;
        virtual ~Base() = default;
    };

    class DetectorReadHits : public Base {
    public:
        virtual void ApplyTo(const readhits_t& hits, extrahits_t& extrahits) = 0;
    };

    class Clusters : public Base {
    public:
        virtual void ApplyTo(clusters_t& clusters) = 0;
    };
};


/**
 * @brief The Updateable_traits class
 *
 * used by the UpdateableManager
 */
class Updateable_traits {
public:
    virtual std::vector<std::list<TID>> GetChangePoints() const = 0;
    virtual void Update(std::size_t index, const TID& id) = 0;
};

} // namespace ant
