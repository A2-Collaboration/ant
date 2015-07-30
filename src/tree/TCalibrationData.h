#ifndef ANT_TCALIBRATIONDATA_H
#define ANT_TCALIBRATIONDATA_H

#include "TDataRecord.h"

#include <string>
#include <ctime>

#define ANT_CALIBRATION_DATA_VERSION 1

namespace ant {

#ifndef __CINT__
struct TCalibrationData : printable_traits
        #else
struct TCalibrationData
        #endif
{
    std::string Author;
    std::string Comment;

    std::int64_t TimeStamp;

    std::string CalibrationID;

    TID FirstID;
    TID LastID;

    typedef TKeyValue<double> TCalibrationEntry;
    std::vector<TCalibrationEntry> Data;

    TCalibrationData() :
        Author(),
        Comment(),
        TimeStamp(),
        CalibrationID(),
        FirstID(),
        LastID(),
        Data()
    {}

    TCalibrationData(const std::string& calibrationID,const TID& first_id, const TID& last_id) :
        Author(),
        Comment(),
        TimeStamp(),
        CalibrationID(calibrationID),
        FirstID(first_id),
        LastID(last_id),
        Data()
    {}

    //Constructors for readout from trees
#ifndef __CINT__
    TCalibrationData(const std::string& author, const std::string& comment,
                     const std::time_t& time,
                     const std::string& calibrationID, const TID& first_id,
                     const TID& last_id,
                     const std::vector<TCalibrationEntry>& data) :
        Author   (author),
        Comment  (comment),
        TimeStamp(time),
        CalibrationID  (calibrationID),
        FirstID  (first_id),
        LastID   (last_id),
        Data     (data)
    {}
#endif

    virtual ~TCalibrationData() {}


#ifndef __CINT__
    virtual std::ostream& Print( std::ostream& s) const override {
        s << "TCalibrationData generated at " << std::asctime(std::localtime(&TimeStamp)) << std::endl
          << "  CalibrationID:  " << CalibrationID << std::endl
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
