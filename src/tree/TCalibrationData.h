#ifndef ANT_TCALIBRATIONDATA_H
#define ANT_TCALIBRATIONDATA_H

#include "TDataRecord.h"
#include <string>

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

    TCalibrationData() : SetupID(), FirstID(), LastID(), Data() {}

    TCalibrationData( const std::string& setupID,const TID& first_id, const TID& last_id) :
        SetupID(setupID),
        FirstID(first_id),
        LastID(last_id),
        Data()
    {}

#ifndef __CINT__
    TCalibrationData( const std::string& setupID, const TID& first_id, const TID& last_id, const TCalibrationEntry& data) :
        SetupID(setupID),
        FirstID(first_id),
        LastID(last_id),
        Data{data}
    {}
    TCalibrationData( const std::string& setupID, const TID& first_id, const TID& last_id, const std::vector<TCalibrationEntry>& data) :
        SetupID(setupID),
        FirstID(first_id),
        LastID(last_id),
        Data(data)
    {}
#endif

    virtual ~TCalibrationData() {}

    std::string SetupID;

    TID FirstID;
    TID LastID;

    std::vector<TCalibrationEntry> Data;

#ifndef __CINT__
    virtual std::ostream& Print( std::ostream& s) const override {
        s << "TCalibrationData" << std::endl
          << "  SetupID:        " << SetupID << std::endl
          << "  Valid for IDs:  [" << FirstID << ", " << LastID << "]" << std::endl
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
