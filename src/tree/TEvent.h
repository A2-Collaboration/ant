#ifndef ANT_TEVENT_H
#define ANT_TEVENT_H

#include "TDataRecord.h"
#include "TTrack.h"
#include "TCluster.h"
#include "TTaggerHit.h"

#ifndef __CINT__
#include <iomanip>
#include <sstream>
#endif

namespace ant {

struct TEvent : TDataRecord
{
  std::vector<ant::TTrack> Tracks;
  std::vector<ant::TTaggerHit> TaggerHits;

#ifndef __CINT__
  TEvent(const TID& id) :
    TDataRecord(id)
  {}

  virtual std::ostream& Print( std::ostream& s) const override {
    s << "TEvent:\n";
    s << " " << TaggerHits.size() << " Taggerhits:\n";
    for(auto& th : TaggerHits) {
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

  TEvent() : TDataRecord() {}
  virtual ~TEvent() {}
  ClassDef(TEvent, ANT_UNPACKER_ROOT_VERSION)

};

}

#endif // ANT_TEVENT_H
