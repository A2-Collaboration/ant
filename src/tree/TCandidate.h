#pragma once

#include "TDataRecord.h"
#include "TCluster.h"
#include <TVector3.h>
#include <vector>

#ifndef __CINT__
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
  double Energy;
  double Time;
  double Theta;
  double Phi;
  double VetoEnergy;
  double TrackerEnergy;

  std::vector<TCluster> Clusters;

#ifndef __CINT__
  TCandidate(
      double E,
      double t,
      double theta,
      double phi,
      const std::vector<TCluster>& clusters = {},
      double vetoE = 0.0,
      double trackerE = 0.0
      ) :
    Energy(E),
    Time(t),
    Theta(theta),
    Phi(phi),
    VetoEnergy(vetoE),
    TrackerEnergy(trackerE),
    Clusters(clusters)
  {}


  virtual std::ostream& Print( std::ostream& s) const override {
    return s << "TCandidate: " << Clusters.size() << " clusters, Energy=" << Energy << " Theta=" << Theta <<", Phi=" << Phi << " VetoEnergy=" << VetoEnergy << " TrackerEnergy=" << TrackerEnergy;
  }
#endif

  TCandidate() {}
  virtual ~TCandidate() {}
  ClassDef(TCandidate, ANT_UNPACKER_ROOT_VERSION)
};

}
