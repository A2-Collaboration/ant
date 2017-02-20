#include "TEventData.h"

using namespace std;
using namespace ant;

TEventData::TEventData(const TID& id) : ID(id) {}
TEventData::TEventData() {}
TEventData::~TEventData() {}

ostream& TEventData::Print(ostream& s) const {
    s << "ID=" << ID << endl;

    s << ">> DetectorReadHits: n=" << DetectorReadHits.size() << endl;
    for(auto& i: DetectorReadHits)
        s << i << endl;

    s << ">> SlowControls: n=" << SlowControls.size() << endl;
    for(auto& i : SlowControls)
        s << i << endl;

    s << ">> UnpackerMessages: n=" << UnpackerMessages.size() << endl;
    for(auto& i : UnpackerMessages)
        s << i << endl;

    s << ">> TaggerHits: n=" << TaggerHits.size() << endl;
    for(auto& i : TaggerHits)
        s << i << endl;

    s << ">> " << Trigger << endl;
    s << ">> " << Target << endl;

    s << ">> Clusters: n=" << Clusters.size() << endl;
    for(auto& i : Clusters)
        s << i; // clusters are multiline

    s << ">> Candidates: n=" << Candidates.size() << endl;
    for(auto& i : Candidates)
        s << i << endl;

    return s;
}

void TEventData::ClearDetectorReadHits()
{
    DetectorReadHits.resize(0);
}