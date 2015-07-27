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
class CandidateBuilder;
class AdaptorTClusterHit;
class Clustering;
class UpdateableManager;
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

    void ApplyHooksToReadHits(
            TDetectorRead& detectorRead,
            sorted_bydetectortype_t<TDetectorReadHit*>& sorted_readhits);

    void BuildHits(
            sorted_bydetectortype_t<TDetectorReadHit*>&& sorted_readhits,
            sorted_bydetectortype_t<reconstruct::AdaptorTClusterHit>& sorted_clusterhits,
            TTagger& event_tagger
            );

    void HandleTagger(
            const std::shared_ptr<TaggerDetector_t>& taggerdetector,
            std::list<TDetectorReadHit*> readhits,
            TTagger& event_tagger);

    void BuildClusters(
            sorted_bydetectortype_t<reconstruct::AdaptorTClusterHit>&& sorted_clusterhits,
            sorted_bydetectortype_t<TCluster>& sorted_clusters,
                       std::vector<TCluster>& insane_clusters
            );

    template<typename T>
    using shared_ptr_list = std::list< std::shared_ptr<T> >;

    using sorted_detectors_t = std::map<Detector_t::Type_t,  std::shared_ptr<Detector_t> >;
    sorted_detectors_t sorted_detectors;

    shared_ptr_list<ReconstructHook::DetectorReadHits> hooks_readhits;
    shared_ptr_list<ReconstructHook::Clusters>         hooks_clusters;


    std::unique_ptr<reconstruct::CandidateBuilder>  candidatebuilder;
    std::unique_ptr<reconstruct::Clustering>        clustering;
    std::unique_ptr<reconstruct::UpdateableManager> updateablemanager;

    template<typename T, typename Base>
    void AddToSharedPtrList(const std::shared_ptr<Base> base, shared_ptr_list<T>& list) {
        const auto& ptr = std::dynamic_pointer_cast<T, Base>(base);
        if(ptr != nullptr)
            list.emplace_back(std::move(ptr));
    }
};

}
