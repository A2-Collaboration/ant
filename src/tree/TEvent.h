#pragma once

#include "TDataRecord.h"
#include "TCandidate.h"
#include "TCluster.h"
#include "TTagger.h"

#ifndef __CINT__
#include <iomanip>
#include <sstream>
#endif

namespace ant {

struct TEvent : TDataRecord
{
    typedef std::vector<ant::TCandidate> candidates_t;
    std::vector<ant::TCandidate> Candidates;
    TTagger Tagger;
    typedef std::vector<ant::TCluster> clusters_t;
    std::vector<ant::TCluster> AllClusters;

    void Clear() {
        Candidates.resize(0);
        Tagger.Clear();
        AllClusters.resize(0);
    }

#ifndef __CINT__
    TEvent(const TID& id) :
        TDataRecord(id)
    {}

    virtual std::ostream& Print( std::ostream& s) const override {
        s << "TEvent ID=" << ID << '\n';
        s << " " << Tagger.Hits.size() << " Taggerhits:\n";
        for(auto& th : Tagger.Hits) {
            s << "  " << th << "\n";
        }
        s << " " << Candidates.size() << " Candidates:\n";
        for(auto& cand : Candidates) {
            s << "  " << cand << "\n";
            for(auto& c : cand.Clusters) {
                s << "   " << c << "\n";
                for(auto& h : c.Hits) {
                    s << "    " << h << "\n";
                }
            }
        }
        return s;
    }
#endif

    TEvent() : TDataRecord(), Candidates(), Tagger() {}
    virtual ~TEvent() {}
    ClassDef(TEvent, ANT_UNPACKER_ROOT_VERSION)

};

}
