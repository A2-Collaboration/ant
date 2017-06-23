#pragma once

#include <memory>
#include <list>

#include "Reconstruct_traits.h"

namespace ant {

struct TTaggerHit;

namespace reconstruct {
class UpdateableManager;
}

class Reconstruct : public Reconstruct_traits {
public:

    using clustering_t       = std::unique_ptr<const reconstruct::Clustering_traits>;
    using candidatebuilder_t = std::unique_ptr<const reconstruct::CandidateBuilder_traits>;

    static clustering_t       GetDefaultClustering();
    static candidatebuilder_t GetDefaultCandidateBuilder();

    Reconstruct(clustering_t clustering_ = GetDefaultClustering(),
                candidatebuilder_t candidatebuilder_ = GetDefaultCandidateBuilder());

    // this method converts a TDetectorRead
    // into a calibrated TEvent
    virtual void DoReconstruct(TEventData& reconstructed) const override;

    virtual ~Reconstruct();

    class Exception : public std::runtime_error {
        using std::runtime_error::runtime_error; // use base class constructor
    };

protected:

    const bool includeIgnoredElements = false;

    // sorted_readhits is mutable in order to
    using sorted_readhits_t = ReconstructHook::Base::readhits_t;
    mutable sorted_readhits_t sorted_readhits;

    void ApplyHooksToReadHits(std::vector<TDetectorReadHit>& detectorReadHits) const;

    template<typename T>
    using sorted_bydetectortype_t = std::map<Detector_t::Type_t, std::vector< T > >;

    void BuildHits(sorted_bydetectortype_t<TClusterHit>& sorted_clusterhits,
            std::vector<TTaggerHit>& taggerhits
            ) const;

    void HandleTagger(const std::shared_ptr<TaggerDetector_t>& taggerdetector,
            const std::vector<std::reference_wrapper<TDetectorReadHit>>& readhits,
            std::vector<TTaggerHit>& taggerhits) const;

    using sorted_clusterhits_t = ReconstructHook::Base::clusterhits_t;
    using sorted_clusters_t = ReconstructHook::Base::clusters_t;
    void BuildClusters(const sorted_clusterhits_t& sorted_clusterhits,
                       sorted_clusters_t& sorted_clusters) const;

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
        // helper method for construction
        static sorted_detectors_t Build();
    };
    const sorted_detectors_t sorted_detectors;

    template<typename T>
    using shared_ptr_list = std::list< std::shared_ptr<T> >;

    const shared_ptr_list<ReconstructHook::DetectorReadHits> hooks_readhits;
    const shared_ptr_list<ReconstructHook::ClusterHits>      hooks_clusterhits;
    const shared_ptr_list<ReconstructHook::Clusters>         hooks_clusters;
    const shared_ptr_list<ReconstructHook::EventData>        hooks_eventdata;

    const clustering_t       clustering;
    const candidatebuilder_t candidatebuilder;
    const std::unique_ptr<reconstruct::UpdateableManager> updateablemanager;
};

}
