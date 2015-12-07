#include "TaggerHit.h"

using namespace ant::analysis::data;
using namespace std;

std::ostream &TaggerHit::Print(std::ostream &stream) const
{
    stream << "TaggerHit: Channel=" << Channel << " PhotonEnergy=" << PhotonEnergy << " MeV  Time=" << Time << "ns";
    return stream;
}
