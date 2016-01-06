#pragma once

#include "tree/TCandidate.h"

//#include "base/printable.h"
//#include "base/types.h"
//#include "base/Detector_t.h"
#include "analysis/data/Cluster.h"
//#include "TVector3.h"

//#include <ostream>
#include <memory>
//#include <vector>

namespace ant {
namespace analysis {
namespace data {

///**
// * @brief The Candidate class
// */
//struct Candidate : printable_traits
//{

//    mev_t ClusterEnergy;
//    radian_t Theta;
//    radian_t Phi;
//    ns_t Time;
//    clustersize_t ClusterSize;
//    Detector_t::Any_t Detector;
//    mev_t VetoEnergy;
//    mev_t TrackerEnergy;

//    ClusterList Clusters;

//    Candidate(const mev_t& clusterEnergy,
//          const radian_t& theta,
//          const radian_t& phi,
//          const ns_t& time,
//          const clustersize_t& clusterSize,
//          const Detector_t::Any_t& detector,
//          const mev_t& vetoEnergy,
//          const mev_t& trackerEnergy
//          ) :
//        ClusterEnergy(clusterEnergy),
//        Theta(theta),
//        Phi(phi),
//        Time(time),
//        ClusterSize(clusterSize),
//        Detector(detector),
//        VetoEnergy(vetoEnergy),
//        TrackerEnergy(trackerEnergy)
//    {}

//    operator TVector3() const { TVector3 p; p.SetMagThetaPhi(1.0, Theta, Phi); return p; }

//    const Cluster* FindFirstCluster(Detector_t::Any_t detector) const {
//        for(const auto& cl : Clusters) {
//            if(cl.Detector & detector) {
//                return std::addressof(cl);
//            }
//        }
//        return nullptr;
//    }

//    const Cluster* FindCaloCluster() const {
//        for(const auto& cl : Clusters) {
//            if(cl.Detector & (Detector_t::Type_t::CB | Detector_t::Type_t::TAPS)) {
//                return std::addressof(cl);
//            }
//        }
//        return nullptr;
//    }

//    const Cluster* FindVetoCluster() const {
//        if(VetoEnergy > 0.0) {
//            for(const auto& cl : Clusters) {
//                if(cl.Detector & (Detector_t::Type_t::PID | Detector_t::Type_t::TAPSVeto)) {
//                    return std::addressof(cl);
//                }
//            }
//        }
//        return nullptr;
//    }

//    virtual ~Candidate() = default;
//    virtual std::ostream& Print(std::ostream &stream) const override;

//};

using CandidatePtr  = std::shared_ptr<TCandidate>;
using CandidateList = std::vector<CandidatePtr>;

}
}
}
