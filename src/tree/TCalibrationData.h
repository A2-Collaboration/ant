#ifndef ANT_TCALIBRATIONDATA_H
#define ANT_TCALIBRATIONDATA_H

#include "base/printable.h"
#include "TDataRecord.h"

#include "Rtypes.h"

#include <map>


#ifndef __CINT__
#include <iomanip>
#include <sstream>
#endif

#define ANT_CALIBRATION_DATA_VERSION 1

namespace ant {

struct TCalibrationEntry
{
  unsigned Key;
  double   Value;
  virtual ~TCalibrationEntry(){}
  ClassDef(TCalibrationEntry, ANT_CALIBRATION_DATA_VERSION)
};

#ifndef __CINT__
struct TCalibrationData : printable_traits
#else
struct TCalibrationData
#endif
{

  TCalibrationData() {}

  TCalibrationData(const TID& first_id, const TID& last_id) :
    FirstID(first_id),
    LastID(last_id)
  {}

  virtual ~TCalibrationData() {}

  TID FirstID;
  TID LastID;

  std::vector<TCalibrationEntry> Data;

#ifndef __CINT__
  virtual std::ostream& Print( std::ostream& s) const override {
    s << "TCalibrationData:\n";
    s << "  [" << FirstID << ", " << LastID << "]\n";
    return s;
  }
#endif

  ClassDef(TCalibrationData, ANT_CALIBRATION_DATA_VERSION)

};

}

#endif // ANT_TCALIBRATIONDATA_H
