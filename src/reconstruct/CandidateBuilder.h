#pragma once

#include "tree/TCandidate.h"
#include "tree/TCluster.h"

#include "expconfig/ExpConfig.h"

#include <map>
#include <list>
#include <memory>

namespace ant {

namespace expconfig {
namespace detector {
struct CB;
struct PID;
struct TAPS;
struct TAPSVeto;
}
}


namespace reconstruct {

class CandidateBuilder {
public:

    using sorted_clusters_t = std::map<Detector_t::Type_t, TClusterList >;
    using candidates_t = TCandidateList;
    using clusters_t = TClusterList;

protected:

    /**
     * @brief option_allowSingleVetoClusters: Make unmached Veto (PID/TAPSVeto) clusters into individual candidates.
     *        To detect partiles that get stuck in the plastic scintillator.
     * @todo Make this configurable via ant::expconfig.
     */
    bool option_allowSingleVetoClusters = false;

    const std::shared_ptr<const expconfig::detector::CB>  cb;
    const std::shared_ptr<const expconfig::detector::PID> pid;
    const std::shared_ptr<const expconfig::detector::TAPS> taps;
    const std::shared_ptr<const expconfig::detector::TAPSVeto> tapsveto;

    const ExpConfig::Setup::candidatebuilder_config_t config;

    void Build_PID_CB(
            sorted_clusters_t& sorted_clusters,
            candidates_t& candidates, clusters_t& all_clusters
            ) const;

    void Build_TAPS_Veto(
            sorted_clusters_t& sorted_clusters,
            candidates_t& candidates, clusters_t& all_clusters
            ) const;

    void Catchall(
            sorted_clusters_t& sorted_clusters,
            candidates_t& candidates, clusters_t& all_clusters
            ) const;

    virtual void BuildCandidates(
            sorted_clusters_t& sorted_clusters,
            candidates_t& candidates,  clusters_t& all_clusters) const;

public:

    CandidateBuilder(const std::shared_ptr<ExpConfig::Setup>& setup);
    virtual ~CandidateBuilder() = default;

    // this method shall fill the TEvent reference
    // with tracks built from the given sorted clusters
    /// \todo make this method abstract and create proper derived Candidate builders
    virtual void Build(
            sorted_clusters_t sorted_clusters,
            candidates_t& candidates,
            clusters_t& all_clusters
            ) const;
};

}} // namespace ant::reconstruct
