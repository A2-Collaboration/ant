#include "TParticle.h"
#include "base/std_ext/math.h"

using namespace ant;
using namespace ant::std_ext;

TParticle::TParticle(const ParticleTypeDatabase::Type& type_, mev_t Ek_, radian_t theta_, radian_t phi_) :
  type(std::addressof(type_))
{
    const mev_t E = Ek_ + type->Mass();
    const mev_t p = sqrt( sqr(E) - sqr(type->Mass()) );

    this->SetPxPyPzE(
                p*TMath::Sin(theta_)*TMath::Cos(phi_),
                p*TMath::Sin(theta_)*TMath::Sin(phi_),
                p*TMath::Cos(theta_),
                E
                );
}

void TParticle::ChangeType(const ParticleTypeDatabase::Type& newtype)
{
    // recalculate Lorentz vector
    const mev_t newE = Ek() + newtype.Mass();
    const mev_t newP = sqrt( sqr(newE) - sqr(newtype.Mass()) );

    SetRho(newP);
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
        stream << " Candidate=" << *Candidate;
    }
    return stream;
}
