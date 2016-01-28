#pragma once

#include "TCluster.h"

#include "base/Detector_t.h"

#include <TVector3.h>

#include <vector>
#include <iomanip>
#include <sstream>


namespace ant {

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

    std::vector<TCluster> Clusters;

    TCandidate(
            Detector_t::Any_t detector,
            double caloE,
            double theta,
            double phi,
            double time,
            unsigned clusterSize,
            double vetoE,
            double trackerE,
            const std::vector<TCluster>& clusters
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

    operator TVector3() const { TVector3 p; p.SetMagThetaPhi(1.0, Theta, Phi); return p; }

    const TCluster* FindFirstCluster(Detector_t::Any_t detector) const {
        for(const auto& cl : Clusters) {
            if(cl.GetDetectorType() & detector) {
                return std::addressof(cl);
            }
        }
        return nullptr;
    }

    const TCluster* FindCaloCluster() const {
        return FindFirstCluster(Detector_t::Type_t::CB | Detector_t::Type_t::TAPS);
    }

    const TCluster* FindVetoCluster() const {
        if(VetoEnergy == 0.0)
            return nullptr;
        return FindFirstCluster(Detector_t::Type_t::PID | Detector_t::Type_t::TAPSVeto);
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
