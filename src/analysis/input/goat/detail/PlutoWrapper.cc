#include "PlutoWrapper.h"

std::ostream& operator <<(std::ostream &stream, const PParticle *p)
{
    stream << "PParticle Type=" << GetParticleName(p) << "(" << p->ID() << ") ParentIndex=" << p->GetParentId()
           << " ParentType=(" << p->GetParentId() << ")";
    return stream;
}


ostream& operator <<(ostream &stream, const PParticle &p)
{
    stream << &p;
    return stream;
}


string GetParticleName(const PParticle *p)
{
    // Unfortunately PParticle::Name() to get the string repesentation of the particle type in
    // pluto is not const. Working around this... (urhg).
    PParticle* ncp = const_cast<PParticle*>(p);
    return ncp->Name();
}


void PrintParticleTable(ostream & stream, const std::vector<PParticle *> plist)
{
    for(std::size_t i=0; i<plist.size(); ++i) {
        stream << i << " " << plist.at(i) << "\n";
    }
}
