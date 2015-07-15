#ifndef ANT_TEVENT_H
#define ANT_TEVENT_H

#include "TDataRecord.h"

#ifndef __CINT__
#include <iomanip>
#include <sstream>
#endif

namespace ant {

struct TEvent : TDataRecord
{
  TEvent(const TDataRecord::ID_t& id) :
    TDataRecord(id)
  {}



#ifndef __CINT__
  virtual std::ostream& Print( std::ostream& s) const override {
    return s << "TEvent";
  }
#endif

  TEvent() : TDataRecord() {}
  ClassDef(TEvent, ANT_UNPACKER_ROOT_VERSION)

};

}

#endif // ANT_TEVENT_H
