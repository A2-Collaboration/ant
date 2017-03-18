#include "TEventData.h"

using namespace std;
using namespace ant;

TEventData::TEventData(const TID& id) : ID(id) {}
TEventData::TEventData() = default;

namespace ant {
ostream& operator<<(ostream& s, const TEventData& o) {
    s << "ID=" << o.ID << endl;

    s << ">> DetectorReadHits: n=" << o.DetectorReadHits.size() << endl;
    for(auto& i: o.DetectorReadHits)
        s << i << endl;

    s << ">> SlowControls: n=" << o.SlowControls.size() << endl;
    for(auto& i : o.SlowControls)
        s << i << endl;

    s << ">> UnpackerMessages: n=" << o.UnpackerMessages.size() << endl;
    for(auto& i : o.UnpackerMessages)
        s << i << endl;

    s << ">> TaggerHits: n=" << o.TaggerHits.size() << endl;
    for(auto& i : o.TaggerHits)
        s << i << endl;

    s << ">> " << o.Trigger << endl;
    s << ">> " << o.Target << endl;

    s << ">> Clusters: n=" << o.Clusters.size() << endl;
    for(auto& i : o.Clusters)
        s << i; // clusters are multiline

    s << ">> Candidates: n=" << o.Candidates.size() << endl;
    for(auto& i : o.Candidates)
        s << i << endl;

    return s;
}
} // namespace ant

void TEventData::ClearDetectorReadHits()
{
    DetectorReadHits.resize(0);
}