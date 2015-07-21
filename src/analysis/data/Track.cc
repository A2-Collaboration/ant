#include "Track.h"

using namespace std;
using namespace ant;

ostream& Track::Print(ostream &stream) const
{
    stream << "Track "
           << " ClusterEnergy=" << ClusterEnergy()
           << " Theta=" << Theta()
           << " Phi=" << Phi()
           << " Time=" << Time()
           << " ClusterSize=" << ClusterSize()
           << " Detectors=" << Detector()
           << " VetoEnergy=" << VetoEnergy()
           << " TrackerEnergy=" << TrackerEnergy();
    return stream;
}
