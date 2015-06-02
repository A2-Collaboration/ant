#include "Particle.h"
#include "TMath.h"

template <class T>
T square(const T& a) { return a*a; }

using namespace ant;

Particle::Particle(const ParticleTypeDatabase::Type &_type, mev_t _Ek, radian_t _theta, radian_t _phi):
    type(&_type)
{
    const mev_t E = _Ek + type->Mass();
    const mev_t p = sqrt( square(E) - square(type->Mass()) );

    //TODO: fix. This might be inefficeint...

    TVector3 pv(1,0,0);

    pv.SetMagThetaPhi(p,_theta,_phi);

    SetLorentzVector(TLorentzVector(pv, E));
}

void Particle::ChangeType(const ParticleTypeDatabase::Type &newtype)
{
    // recalculate Lorentz vector
    const mev_t newE = Ek() + newtype.Mass();
    const mev_t newP = sqrt( square(newE) - square(newtype.Mass()) );

    TVector3 pv = Vect();
    pv.SetMag(newP);

    SetVect(pv);
    SetE(newE);
    type = &newtype;

}

std::ostream &Particle::Print(std::ostream &stream) const
{
    stream << "Particle " << Type().Name();
    stream << " IM=" << M();
    stream << " E=" << E();
    stream << " Theta=" << Theta();
    stream << " Phi=" << Phi();
    for( auto& track : tracks ) {
        stream << "\t" << track << "\n";
    }
    return stream;
}
