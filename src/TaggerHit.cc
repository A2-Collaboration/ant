#include "TaggerHit.h"

using namespace ant;
using namespace std;

std::ostream &TaggerHit::Print(std::ostream &stream) const
{
    stream << "TaggerHit: Channel=" << channel << " PhotonEnergy=" << photon_energy << " MeV  Time=" << time << "ns";
    return stream;
}
