#include "TParticle.h"
#include "base/std_ext/math.h"

using namespace ant;
using namespace ant::std_ext;

TParticle::TParticle(const ParticleTypeDatabase::Type& _type, mev_t _Ek, radian_t _theta, radian_t _phi) :
  type(std::addressof(_type))
{
    const mev_t E = _Ek + type->Mass();
    const mev_t p = sqrt( sqr(E) - sqr(type->Mass()) );

    /// \bug This might be inefficient...

    TVector3 pv(1,0,0);

    pv.SetMagThetaPhi(p,_theta,_phi);

    *(reinterpret_cast<TLorentzVector*>(this)) = TLorentzVector(pv, E);
}

void TParticle::ChangeType(const ParticleTypeDatabase::Type& newtype)
{
    // recalculate Lorentz vector
    const mev_t newE = Ek() + newtype.Mass();
    const mev_t newP = sqrt( sqr(newE) - sqr(newtype.Mass()) );

    TVector3 pv = Vect();
    pv.SetMag(newP);

    SetVect(pv);
    SetE(newE);
    type = &newtype;

}

std::ostream &TParticle::Print(std::ostream &stream) const
{
    stream << "TParticle " << Type().Name();
    stream << " IM=" << M();
    stream << " E=" << E();
    stream << " Theta=" << Theta();
    stream << " Phi=" << Phi();
    if(Candidate) {
        stream << " Candidate=" << *Candidate << "\n";
    }
    return stream;
}
