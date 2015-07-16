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
  TEvent() : TDataRecord() {}
  TEvent(const TID& id) :
    TDataRecord(id)
  {}
  virtual ~TEvent() {}

  std::vector<ant::TTrack> Tracks;
  std::vector<ant::TTaggerHit> TaggerHits;

#ifndef __CINT__
  virtual std::ostream& Print( std::ostream& s) const override {
    s << "TEvent:";
    s << " " << TaggerHits.size() << " Taggerhits";
    s << " " << Tracks.size() << " Tracks";
    return s;
  }
#endif

  ClassDef(TEvent, ANT_UNPACKER_ROOT_VERSION)

};

}

#endif // ANT_TEVENT_H
