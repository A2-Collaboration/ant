#pragma once

#include <memory>
#include <list>

#include "Reconstruct_traits.h"

namespace ant {

struct TEvent;
struct THeaderInfo;
struct TDetectorRead;
struct TCluster;
struct TTagger;

namespace reconstruct {
class CandidateBuilder;
struct AdaptorTClusterHit;
class Clustering;
class UpdateableManager;
}

class Reconstruct : public Reconstruct_traits {
    // used in test/reconstruct/TestReconstruct.cc
    // to test the private methods interplay
    friend struct ReconstructTester;

public:
    Reconstruct();

    // this method converts a TDetectorRead
    // into a calibrated TEvent
    virtual void DoReconstruct(TEvent::Data& reconstructed) override;

    ~Reconstruct();

    class Exception : public std::runtime_error {
        using std::runtime_error::runtime_error; // use base class constructor
    };

    template<typename T>
    using sorted_bydetectortype_t = std::map<Detector_t::Type_t, std::list< T > >;

private:

    bool initialized = false;
    virtual void Initialize(const TID& tid);

    using sorted_readhits_t = ReconstructHook::Base::readhits_t;
    sorted_readhits_t sorted_readhits;

    void ApplyHooksToReadHits(std::vector<TDetectorReadHit>& detectorReadHits);

    void BuildHits(
            sorted_bydetectortype_t<reconstruct::AdaptorTClusterHit>& sorted_clusterhits,
            TTagger& event_tagger
            );

    void HandleTagger(const std::shared_ptr<TaggerDetector_t>& taggerdetector,
            const std::vector<TDetectorReadHit*>& readhits,
            TTagger& event_tagger);

    void BuildClusters(sorted_bydetectortype_t<reconstruct::AdaptorTClusterHit>&& sorted_clusterhits,
            sorted_bydetectortype_t<TClusterPtr>& sorted_clusters);



    // little helper class which stores the upcasted versions of shared_ptr
    // to Detector_t instances
    struct detector_ptr_t {
        std::shared_ptr<Detector_t> Detector;
        std::shared_ptr<TaggerDetector_t> TaggerDetector; // might be nullptr
        std::shared_ptr<ClusterDetector_t> ClusterDetector; // might be nullptr

        detector_ptr_t(const std::shared_ptr<Detector_t>& detector) :
            Detector(detector),
            TaggerDetector(std::dynamic_pointer_cast<TaggerDetector_t>(detector)),
            ClusterDetector(std::dynamic_pointer_cast<ClusterDetector_t>(detector))
        {
            if(TaggerDetector != nullptr && ClusterDetector != nullptr) {
                throw Exception("Found detector which is both clustering and tagging, not supported");
            }
        }
        // implicit conversion to simple base class Detector pointer
        operator std::shared_ptr<Detector_t>() const { return Detector; }
    };
    struct sorted_detectors_t : std::map<Detector_t::Type_t,  detector_ptr_t > {
        // use the implicit conversions of detector_ptr_t
        operator std::map<Detector_t::Type_t,  std::shared_ptr<Detector_t> >() const {
            return std::map<Detector_t::Type_t,  std::shared_ptr<Detector_t> >(begin(), end());
        }
    };
    sorted_detectors_t sorted_detectors;

    template<typename T>
    using shared_ptr_list = std::list< std::shared_ptr<T> >;
    shared_ptr_list<ReconstructHook::DetectorReadHits> hooks_readhits;
    shared_ptr_list<ReconstructHook::ClusterHits>      hooks_clusterhits;
    shared_ptr_list<ReconstructHook::Clusters>         hooks_clusters;
    shared_ptr_list<ReconstructHook::EventData>        hooks_eventdata;


    std::unique_ptr<reconstruct::CandidateBuilder>  candidatebuilder;
    std::unique_ptr<reconstruct::Clustering>        clustering;
    std::unique_ptr<reconstruct::UpdateableManager> updateablemanager;
};

}
