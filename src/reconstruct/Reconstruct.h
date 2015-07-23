#pragma once

#include <memory>
#include <list>

#include "Reconstruct_traits.h"

namespace ant {

class TEvent;
class THeaderInfo;
class TDetectorRead;
class TCluster;
class TTagger;

namespace reconstruct {
class TrackBuilder;
class HitWithEnergy_t;
class Clustering;
}

class Reconstruct {
    // used in test/reconstruct/TestReconstruct.cc
    // to test the private methods interplay
    friend class ReconstructTester;

public:
    // You can only create the reconstruct machinery
    // if it's able to find its config. For this, it needs the
    // some THeaderInfo object
    Reconstruct(const THeaderInfo& headerInfo);

    // this method converts a TDetectorRead
    // into a calibrated TEvent
    std::unique_ptr<TEvent> DoReconstruct(TDetectorRead& detectorRead);

    ~Reconstruct();

private:

    template<typename T>
    using sorted_bydetectortype_t = std::map<Detector_t::Type_t, std::list< T > >;

    void BuildHits(
            sorted_bydetectortype_t<TDetectorReadHit*>&& sorted_readhits,
            sorted_bydetectortype_t<reconstruct::HitWithEnergy_t>& sorted_clusterhits,
            TTagger& event_tagger
            );
    void BuildClusters(sorted_bydetectortype_t<reconstruct::HitWithEnergy_t>&& sorted_clusterhits,
            sorted_bydetectortype_t<TCluster>& sorted_clusters
            );

    template<typename T>
    using shared_ptr_list = std::list< std::shared_ptr<T> >;

    using sorted_detectors_t = std::map<Detector_t::Type_t, std::shared_ptr<Detector_t> >;


    shared_ptr_list<CalibrationApply_traits>  calibrations;
    shared_ptr_list<Updateable_traits>        updateables;
    sorted_detectors_t                        sorted_detectors;

    std::unique_ptr<reconstruct::TrackBuilder> trackbuilder;
    std::unique_ptr<reconstruct::Clustering>   clustering;


};

}
