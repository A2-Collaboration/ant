#include "TSimpleParticle.h"
#include "base/std_ext/math.h"

using namespace ant;
using namespace ant::std_ext;


TSimpleParticle::TSimpleParticle(const TParticle& particle):
    LorentzVec(particle)
{
    VetoE       = particle.Candidate->VetoEnergy;
    ShortE      = particle.Candidate->FindCaloCluster()->ShortEnergy;
    ClusterSize = particle.Candidate->ClusterSize;
}

namespace ant {
std::ostream& operator<<(std::ostream &stream, const TSimpleParticle& o)
{
    stream << " IM=" << o.M();
    stream << " E=" << o.E;
    stream << " Theta=" << o.Theta();
    stream << " Phi=" << o.Phi();
    stream << " VetoE=" << o.VetoE;
    stream << " ClusterSize=" << o.ClusterSize;
    return stream;
}

}
