#ifndef ANT_TTRACK_H
#define ANT_TTRACK_H

#include "TDataRecord.h"
#include "TCluster.h"
#include <TVector3.h>
#include <vector>

#ifndef __CINT__
#include <iomanip>
#include <sstream>
#include <limits> // for NaN
#endif

namespace ant {

#ifndef __CINT__
struct TTrack: ant::printable_traits
#else
struct TTrack
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
  TTrack(
      double E,
      double t,
      double theta,
      double phi,
      const std::vector<TCluster>& clusters = {},
      double vetoE = std::numeric_limits<double>::quiet_NaN(),
      double trackerE = std::numeric_limits<double>::quiet_NaN()
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
    return s << "TTrack: " << Clusters.size() << " clusters Energy=" << Energy << " Theta=" << Theta <<", Phi=" << Phi << " VetoEnergy=" << VetoEnergy << " TrackerEnergy=" << TrackerEnergy;
  }
#endif

  TTrack() {}
  virtual ~TTrack() {}
  ClassDef(TTrack, ANT_UNPACKER_ROOT_VERSION)
};

}

#endif // ANT_TTRACK_H
