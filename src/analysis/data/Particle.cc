#include "Particle.h"
#include "TMath.h"

template <class T>
T square(const T& a) { return a*a; }

using namespace ant;
using namespace ant::analysis::data;

Particle::Particle(const ParticleTypeDatabase::Type &_type, mev_t _Ek, radian_t _theta, radian_t _phi) :
  type(&_type),
  parents(),
  daughters(),
  candidates()
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
    for( auto& candidate : candidates ) {
        stream << "\t" << *candidate << "\n";
    }
    return stream;
}


void Particle::RecPrint(const ParticlePtr &p, std::ostream &stream)
{

    stream << p->Type().PrintName() << " ";

    if(! p->Daughters().empty()) {

        if(p->Type() == ParticleTypeDatabase::BeamTarget) {
            stream << "#rightarrow ";
            for(auto& d : p->Daughters()) {
                RecPrint(d, stream);
            }
            stream << " ";
        } else {
            stream << "[ ";
            for(auto& d : p->Daughters()) {
                RecPrint(d, stream);
            }
            stream << "] ";
        }

    }
}
