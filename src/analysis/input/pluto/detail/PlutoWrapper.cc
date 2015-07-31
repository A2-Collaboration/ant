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
    return GetParticleName(p->ID());
}

string GetParticleName(const int id)
{
    if(id == 14001) return "gp";
    return makeStaticData()->GetParticleName(id);
}

void PrintParticleTable(ostream & stream, const std::vector<const PParticle *> plist)
{
    for(std::size_t i=0; i<plist.size(); ++i) {
        stream << i << " " << plist.at(i) << "\n";
    }
}
