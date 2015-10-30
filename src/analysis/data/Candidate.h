#pragma once

#include "base/printable.h"
#include "base/types.h"
#include "base/Detector_t.h"
#include "analysis/data/Cluster.h"
#include "TVector3.h"

#include <ostream>
#include <memory>
#include <vector>

namespace ant {
namespace analysis {
namespace data {

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
    Detector_t::Any_t detector;
    mev_t vetoEnergy;
    mev_t trackerEnergy;

    ClusterList Clusters;

    Candidate(const mev_t& _clusterEnergy,
          const radian_t& _theta,
          const radian_t& _phi,
          const ns_t& _time,
          const clustersize_t& _clusterSize,
          const Detector_t::Any_t& _detector,
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
    Detector_t::Any_t Detector() const { return detector; }
    mev_t VetoEnergy() const { return vetoEnergy; }
    mev_t TrackerEnergy() const { return trackerEnergy; }
    operator TVector3 () const   { TVector3 p; p.SetMagThetaPhi(1.0, Theta(), Phi()); return p; }

    const Cluster* FindFirstCluster(Detector_t::Any_t detector) {
        for(const auto& cl : Clusters) {
            if(cl.Detector & detector) {
                return std::addressof(cl);
            }
        }
        return nullptr;
    }

    const Cluster* FindCaloCluster() {
        for(const auto& cl : Clusters) {
            if(cl.Detector & (Detector_t::Type_t::CB | Detector_t::Type_t::TAPS)) {
                return std::addressof(cl);
            }
        }
        return nullptr;
    }

    const Cluster* FindVetoCluster() {
        if(VetoEnergy() > 0.0) {
            for(const auto& cl : Clusters) {
                if(cl.Detector & (Detector_t::Type_t::PID | Detector_t::Type_t::TAPSVeto)) {
                    return std::addressof(cl);
                }
            }
        }
        return nullptr;
    }

    virtual std::ostream& Print(std::ostream &stream) const override;

};

using CandidatePtr  = std::shared_ptr<Candidate>;
using CandidateList = std::vector<CandidatePtr>;

}
}
}
