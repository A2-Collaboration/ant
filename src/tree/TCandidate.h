#pragma once

#include "TCluster.h"

#include "base/Detector_t.h"

#include <TVector3.h>

#include <vector>
#include <iomanip>
#include <sstream>
#include <memory>
#include <algorithm>


namespace ant {

struct TCandidate;

using TCandidatePtr = std::shared_ptr<TCandidate>;
using TCandidateList = std::vector<TCandidatePtr>;

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

    TClusterList Clusters;

    TCandidate(
            Detector_t::Any_t detector,
            double caloE,
            double theta,
            double phi,
            double time,
            unsigned clusterSize,
            double vetoE,
            double trackerE,
            TClusterList clusters
            ) :
        Detector(detector),
        CaloEnergy(caloE),
        Theta(theta),
        Phi(phi),
        Time(time),
        ClusterSize(clusterSize),
        VetoEnergy(vetoE),
        TrackerEnergy(trackerE),
        Clusters(std::move(clusters))
    {}

    template<class Archive>
    void serialize(Archive& archive) {
        archive(Detector, CaloEnergy, Theta, Phi, Time, ClusterSize, VetoEnergy, TrackerEnergy, Clusters);
    }

    operator TVector3() const { TVector3 p; p.SetMagThetaPhi(1.0, Theta, Phi); return p; }

    TClusterList::const_iterator FindFirstCluster(Detector_t::Any_t detector) const {
        return std::find_if(Clusters.begin(), Clusters.end(), [detector] (const TCluster& cl) {
            return cl.DetectorType & detector;
        });
    }

    TClusterList::const_iterator FindCaloCluster() const {
        return FindFirstCluster(Detector_t::Any_t::Calo);
    }

    TClusterList::const_iterator FindVetoCluster() const {
        return FindFirstCluster(Detector_t::Any_t::Veto);
    }

    virtual std::ostream& Print( std::ostream& s) const override {
        return s << "TCandidate: " << Clusters.size() << " clusters, Detector=" << Detector
                 << " CaloEnergy=" << CaloEnergy
                 << " Theta=" << Theta <<", Phi=" << Phi
                 << " ClusterSize=" << ClusterSize << " Time=" << Time
                 << " VetoEnergy=" << VetoEnergy << " TrackerEnergy=" << TrackerEnergy;
    }

    TCandidate() : Detector(Detector_t::Any_t::None) {}
    virtual ~TCandidate() {}
};


}
