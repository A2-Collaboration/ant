#include "Candidate.h"

using namespace std;
using namespace ant::analysis::data;

ostream& Candidate::Print(ostream &stream) const
{
    stream << "Candidate "
           << " ClusterEnergy=" << ClusterEnergy
           << " Theta=" << Theta
           << " Phi=" << Phi
           << " Time=" << Time
           << " ClusterSize=" << ClusterSize
           << " Detectors=" << Detector
           << " VetoEnergy=" << VetoEnergy
           << " TrackerEnergy=" << TrackerEnergy;
    return stream;
}
