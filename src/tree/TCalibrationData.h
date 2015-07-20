#ifndef ANT_TCALIBRATIONDATA_H
#define ANT_TCALIBRATIONDATA_H

#include "TDataRecord.h"

#define ANT_CALIBRATION_DATA_VERSION 1

namespace ant {

struct TCalibrationEntry
{
  unsigned Key;
  double   Value;

  TCalibrationEntry() : Key(), Value() {}
  TCalibrationEntry(unsigned key, double value) : Key(key), Value(value) {}

  virtual ~TCalibrationEntry(){}
  ClassDef(TCalibrationEntry, ANT_CALIBRATION_DATA_VERSION)
};

#ifndef __CINT__
struct TCalibrationData : printable_traits
#else
struct TCalibrationData
#endif
{

  TCalibrationData() : FirstID(), LastID(), Data() {}

  TCalibrationData(const TID& first_id, const TID& last_id) :
    FirstID(first_id),
    LastID(last_id),
    Data()
  {}

  #ifndef __CINT__
  TCalibrationData(const TID& first_id, const TID& last_id, const std::vector<TCalibrationEntry>& data) :
    FirstID(first_id),
    LastID(last_id),
    Data(data)
  {}
  #endif

  virtual ~TCalibrationData() {}

  TID FirstID;
  TID LastID;

  std::vector<TCalibrationEntry> Data;

#ifndef __CINT__
  virtual std::ostream& Print( std::ostream& s) const override {
    s << "TCalibrationData:\n";
    s << "  Valid for IDs:  [" << FirstID << ", " << LastID << "]" << std::endl
      << "  Data:" << std::endl;
    for (auto& entry: Data){
        s << "  " << entry.Key << "  " << entry.Value << std::endl;
    }
    return s;
  }
#endif

  ClassDef(TCalibrationData, ANT_CALIBRATION_DATA_VERSION)

};

}

#endif // ANT_TCALIBRATIONDATA_H
