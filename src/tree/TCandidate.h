#pragma once

#include "TCluster.h"

#include "base/Detector_t.h"

#include <TVector3.h>

#include <vector>
#include <iomanip>
#include <sstream>
#include <memory>


namespace ant {

struct TCandidate;

using TCandidatePtr = std::shared_ptr<TCandidate>;

struct TCandidate : printable_traits
{
    Detector_t::Any_t Detector;
    double CaloEnergy;
    double Theta;
    double Phi;
    double Time;
    std::uint16_t ClusterSize;

    double VetoEnergy;
    double TrackerEnergy;

    std::vector<TClusterPtr> Clusters;

    TCandidate(
            Detector_t::Any_t detector,
            double caloE,
            double theta,
            double phi,
            double time,
            unsigned clusterSize,
            double vetoE,
            double trackerE,
            const std::vector<TClusterPtr>& clusters
            ) :
        Detector(detector),
        CaloEnergy(caloE),
        Theta(theta),
        Phi(phi),
        Time(time),
        ClusterSize(clusterSize),
        VetoEnergy(vetoE),
        TrackerEnergy(trackerE),
        Clusters(clusters)
    {}

    template<class Archive>
    void serialize(Archive& archive) {
        archive(Detector, CaloEnergy, Theta, Phi, Time, ClusterSize, VetoEnergy, TrackerEnergy, Clusters);
    }

    operator TVector3() const { TVector3 p; p.SetMagThetaPhi(1.0, Theta, Phi); return p; }

    TClusterPtr FindFirstCluster(Detector_t::Any_t detector) const {
        for(const auto& cl : Clusters) {
            if(cl->DetectorType & detector) {
                return cl;
            }
        }
        return nullptr;
    }

    TClusterPtr FindCaloCluster() const {
        return FindFirstCluster(Detector_t::Type_t::CB | Detector_t::Type_t::TAPS);
    }

    TClusterPtr FindVetoCluster() const {
        return FindFirstCluster(Detector_t::Any_t::Veto);
    }

    virtual std::ostream& Print( std::ostream& s) const override {
        return s << "TCandidate: " << Clusters.size() << " clusters, CaloEnergy=" << CaloEnergy
                 << " Theta=" << Theta <<", Phi=" << Phi
                 << " ClusterSize=" << ClusterSize << " Time=" << Time
                 << " VetoEnergy=" << VetoEnergy << " TrackerEnergy=" << TrackerEnergy;
    }

    TCandidate() : Detector(Detector_t::Any_t::None) {}
    virtual ~TCandidate() {}
};


}
