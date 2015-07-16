#ifndef ANT_TTRACK_H
#define ANT_TTRACK_H

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
struct TTrack: ant::printable_traits
#else
struct TTrack
#endif
{
    std::vector<ant::TCluster> Clusters;
    double Energy;
    double Time;
    double Theta;
    double Phi;
    double VetoEnergy;
    double TrackerEnergy;

    TTrack(): Energy(0.0), Theta(0.0), Phi(0.0), VetoEnergy(0.0), TrackerEnergy(0.0) {}

    TTrack(double E, double t, double theta, double phi, double vetoE, double trackerE ):
        Energy(E), Time(t), Theta(theta), Phi(phi), VetoEnergy(vetoE), TrackerEnergy(trackerE) {}

    virtual ~TTrack() {}

#ifndef __CINT__
  virtual std::ostream& Print( std::ostream& s) const override {
    return s << "TTrack: " << Clusters.size() << " clusters Theta=" << Theta <<", Phi=" << Phi << " VetoEnergy=" << VetoEnergy << " TrackerEnergy=" << TrackerEnergy;
  }
#endif

    ClassDef(TTrack, ANT_UNPACKER_ROOT_VERSION)

};

}

#endif // ANT_TTRACK_H
