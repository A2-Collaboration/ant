#ifndef ANT_TEVENT_H
#define ANT_TEVENT_H

#include "TDataRecord.h"
#include "TTrack.h"
#include "TCluster.h"
#include "TTagger.h"

#ifndef __CINT__
#include <iomanip>
#include <sstream>
#endif

namespace ant {

struct TEvent : TDataRecord
{
    typedef std::vector<ant::TTrack> tracks_t;
    tracks_t Tracks;
    TTagger Tagger;

#ifndef __CINT__
    TEvent(const TID& id) :
        TDataRecord(id)
    {}

    virtual std::ostream& Print( std::ostream& s) const override {
        s << "TEvent:\n";
        s << " " << Tagger.Hits.size() << " Taggerhits:\n";
        for(auto& th : Tagger.Hits) {
            s << "  " << th << "\n";
        }
        s << " " << Tracks.size() << " Tracks:\n";
        for(auto& t : Tracks) {
            s << "  " << t << "\n";
            for(auto& c : t.Clusters) {
                s << "   " << c << "\n";
                for(auto& h : c.Hits) {
                    s << "    " << h << "\n";
                }
            }
        }
        return s;
    }
#endif

    TEvent() : TDataRecord(), Tracks(), Tagger() {}
    virtual ~TEvent() {}
    ClassDef(TEvent, ANT_UNPACKER_ROOT_VERSION)

};

}

#endif // ANT_TEVENT_H
