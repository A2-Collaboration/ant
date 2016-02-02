#include "TEventData.h"

using namespace std;
using namespace ant;

TEventData::TEventData(const TID& id) : ID(id) {}
TEventData::TEventData() {}
TEventData::~TEventData() {}

const TParticleList TEventData::PTypeList::empty;

ostream& TEventData::Print(ostream& s) const {
    s << "ID=" << ID << endl;

    s << ">> DetectorReadHits" << endl;
    for(auto& i: DetectorReadHits)
        s << i << endl;

    s << ">> SlowControls" << endl;
    for(auto& i : SlowControls)
        s << i << endl;

    s << ">> UnpackerMessages" << endl;
    for(auto& i : UnpackerMessages)
        s << i << endl;

    s << ">> TaggerHits" << endl;
    for(auto& i : TaggerHits)
        s << i << endl;

    s << ">> Clusters" << endl;
    for(auto& i : Clusters)
        s << *i << endl;

    s << ">> Candidates" << endl;
    for(auto& i : Candidates)
        s << *i << endl;

    s << ">> Particles" << endl;
    for(auto& i : Particles.GetAll())
        s << *i << endl;

    return s;
}