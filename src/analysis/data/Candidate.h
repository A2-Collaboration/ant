#pragma once

#include "base/printable.h"
#include "base/types.h"
#include "analysis/Detector.h"
#include "analysis/data/Cluster.h"

#include <ostream>
#include <memory>
#include <vector>

namespace ant {

/**
 * @brief The Track class
 * Representation of GoAT information, with emphasis on
 * physical particle information, not on detector information!
 */
class Candidate: public ant::printable_traits
{
public:
    mev_t clusterEnergy;
    radian_t theta;
    radian_t phi;
    ns_t time;
    clustersize_t clusterSize;
    detector_t detector;
    mev_t vetoEnergy;
    mev_t trackerEnergy;

    std::vector<ant::Cluster> Clusters;

    Candidate(const mev_t& _clusterEnergy,
          const radian_t& _theta,
          const radian_t& _phi,
          const ns_t& _time,
          const clustersize_t& _clusterSize,
          const detector_t& _detector,
          const mev_t& _vetoEnergy,
          const mev_t& _trackerEnergy
          ) :
        clusterEnergy(_clusterEnergy),
        theta(_theta),
        phi(_phi),
        time(_time),
        clusterSize(_clusterSize),
        detector(_detector),
        vetoEnergy(_vetoEnergy),
        trackerEnergy(_trackerEnergy)
    {}


    mev_t ClusterEnergy() const { return clusterEnergy; }
    radian_t Theta() const { return theta; }
    radian_t Phi() const { return phi; }
    ns_t Time() const { return time; }
    clustersize_t ClusterSize() const { return clusterSize; }
    detector_t Detector() const { return detector; }
    mev_t VetoEnergy() const { return vetoEnergy; }
    mev_t TrackerEnergy() const { return trackerEnergy; }

    const Cluster* FindCaloCluster() {
        for(const auto& cl : Clusters) {
            if(cl.Detector & (detector_t::CB | detector_t::TAPS)) {
                return std::addressof(cl);
            }
        }
        return nullptr;
    }

    virtual std::ostream& Print(std::ostream &stream) const override;

};

using CandidatePtr  = std::shared_ptr<Candidate>;
using CandidateList = std::vector<CandidatePtr>;

}
