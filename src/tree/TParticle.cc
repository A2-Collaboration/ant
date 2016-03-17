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

std::ostream &TParticle::Print(std::ostream &stream) const
{
    stream << "TParticle " << Type().Name();
    stream << " IM=" << M();
    stream << " E=" << E;
    stream << " Theta=" << Theta();
    stream << " Phi=" << Phi();
    if(Candidate) {
        stream << " Candidate=" << *Candidate;
    }
    return stream;
}
