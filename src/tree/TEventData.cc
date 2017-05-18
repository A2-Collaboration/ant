#include "TEventData.h"

using namespace std;
using namespace ant;

TEventData::TEventData(const TID& id) : ID(id) {}
TEventData::TEventData() = default;


string _GetDecayString(const TParticleTree_t& particletree)
{
    if(!particletree)
        return "empty_unknown";

    stringstream s;

    // the head is the beam part
    s << particletree->Get()->Type().Name() << " -> ";

    // ignore level==0 since it's the already handled beamparticle
    size_t lastlevel = 1;
    particletree->Map_level([&s, &lastlevel] (const TParticlePtr& p, size_t level) {
        if(level>0) {
            while(lastlevel<level) {
                s << "[ "; lastlevel++;
            }
            while(lastlevel>level) {
                s << "] "; lastlevel--;
            }
            assert(level == lastlevel);
            s << p->Type().Name() << " ";
        }
    });

    while(lastlevel-- > 1)
        s << "] ";

    return s.str();
}

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

    if (o.ParticleTree)
    {
        s << ">> ParticleTree: " << _GetDecayString(o.ParticleTree) << endl;
    }

    return s;
}
} // namespace ant

void TEventData::ClearDetectorReadHits()
{
    DetectorReadHits.resize(0);
}
