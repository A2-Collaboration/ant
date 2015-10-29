#pragma once

#include "base/Detector_t.h"
#include "tree/TCluster.h"
#include "tree/TEvent.h"
#include "expconfig/ExpConfig.h"

#include <map>
#include <list>
#include <memory>

namespace ant {

namespace expconfig {
namespace detector {
class CB;
class PID;
class TAPS;
class TAPSVeto;
}
}


namespace reconstruct {

class CandidateBuilder {
protected:

    /**
     * @brief option_allowSingleVetoClusters: Make unmached Veto (PID/TAPSVeto) clusters into individual candidates.
     *        To detect partiles that get stuck in the plastic scintillator.
     * @todo Make this configurable via ant::expconfig.
     */
    bool option_allowSingleVetoClusters = false;

    std::shared_ptr<expconfig::detector::CB>  cb;
    std::shared_ptr<expconfig::detector::PID> pid;
    std::shared_ptr<expconfig::detector::TAPS> taps;
    std::shared_ptr<expconfig::detector::TAPSVeto> tapsveto;

    const ant::ExpConfig::Reconstruct::candidatebuilder_config_t config;

    void Build_PID_CB(
            std::map<Detector_t::Type_t, std::list< TCluster > >& sorted_clusters,
            TEvent::candidates_t& candidates,
            std::vector<TCluster>& all_clusters
            );

    void Build_TAPS_Veto(
            std::map<Detector_t::Type_t, std::list< TCluster > >& sorted_clusters,
            TEvent::candidates_t& candidates,
            std::vector<TCluster>& all_clusters
            );

    void Catchall(
            std::map<Detector_t::Type_t, std::list< TCluster > >& sorted_clusters,
            TEvent::candidates_t& candidates,
            std::vector<TCluster>& all_clusters
            );

    virtual void BuildCandidates(
            std::map<Detector_t::Type_t, std::list<TCluster> >& sorted_clusters,
            TEvent::candidates_t& candidates,
            std::vector<TCluster>& all_clusters);

public:

    using sorted_detectors_t = std::map<Detector_t::Type_t, std::shared_ptr<Detector_t> >;

    CandidateBuilder(const sorted_detectors_t& sorted_detectors, const std::shared_ptr<ExpConfig::Reconstruct>& _config);
    virtual ~CandidateBuilder() = default;

    // this method shall fill the TEvent reference
    // with tracks built from the given sorted clusters
    /// \todo make this method abstract and create proper derived Candidate builders
    virtual void Build(
            std::map<Detector_t::Type_t, std::list<TCluster> > sorted_clusters,
            TEvent::candidates_t& candidates,
            std::vector<TCluster>& all_clusters
            );
};

}} // namespace ant::reconstruct
