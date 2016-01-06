#pragma once

#include "TDataRecord.h"
#include "TCluster.h"


#include <TVector3.h>

#include <vector>

#ifndef __CINT__
#include "base/Detector_t.h"
#include <iomanip>
#include <sstream>
#endif

namespace ant {

#ifndef __CINT__
struct TCandidate: ant::printable_traits
#else
struct TCandidate
#endif
{
    std::uint64_t Detector;
    double CaloEnergy;
    double Theta;
    double Phi;
    double Time;
    std::uint16_t ClusterSize;

    double VetoEnergy;
    double TrackerEnergy;

    std::vector<TCluster> Clusters;

#ifndef __CINT__
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
        Detector(detector.bitfield),
        CaloEnergy(caloE),
        Theta(theta),
        Phi(phi),
        Time(time),
        ClusterSize(clusterSize),
        VetoEnergy(vetoE),
        TrackerEnergy(trackerE),
        Clusters(clusters)
    {}


    virtual std::ostream& Print( std::ostream& s) const override {
        return s << "TCandidate: " << Clusters.size() << " clusters, CaloEnergy=" << CaloEnergy << " Theta=" << Theta <<", Phi=" << Phi << " VetoEnergy=" << VetoEnergy << " TrackerEnergy=" << TrackerEnergy;
    }
#endif

    TCandidate() {}
    virtual ~TCandidate() {}
    ClassDef(TCandidate, ANT_UNPACKER_ROOT_VERSION)
};

}
