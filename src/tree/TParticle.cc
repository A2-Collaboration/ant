#include "TParticle.h"
#include "base/std_ext/math.h"

using namespace ant;
using namespace ant::std_ext;

TParticle::TParticle(const ParticleTypeDatabase::Type& type_, mev_t Ek_, radian_t theta_, radian_t phi_) :
  LorentzVec(vec3::RThetaPhi(
                 sqrt( sqr(Ek_ + type_.Mass()) - sqr(type_.Mass())), // this could be simplified...
                 theta_,
                 phi_
                 ),
             Ek_ + type_.Mass()),
  type(std::addressof(type_))
{
}

void TParticle::ChangeType(const ParticleTypeDatabase::Type& newtype)
{
    // recalculate Lorentz vector
    E = Ek() + newtype.Mass();
    const mev_t newP = sqrt( sqr(E) - sqr(newtype.Mass()) );
    p.SetR(newP);

    // remember new type
    type = std::addressof(newtype);
}

namespace ant {
std::ostream& operator<<(std::ostream &stream, const TParticle& o)
{
    stream << "TParticle " << o.Type().Name();
    stream << " IM=" << o.M();
    stream << " E=" << o.E;
    stream << " Theta=" << o.Theta();
    stream << " Phi=" << o.Phi();
    if(o.Candidate) {
        stream << " Candidate=" << *o.Candidate;
    }
    return stream;
}
}
