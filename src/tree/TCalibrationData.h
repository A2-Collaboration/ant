#pragma once

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

    bool Extendable; // default is false
    TID FirstID;
    TID LastID;

    typedef TKeyValue<double> TCalibrationEntry;
    std::vector<TCalibrationEntry> Data;

    typedef TKeyValue<std::vector<double>> TFitParameters;
    std::vector<TFitParameters> FitParameters;

#ifndef __CINT__
    TCalibrationData(const std::string& author,
                     const std::string& comment,
                     const std::time_t& time,
                     const std::string& calibrationID,
                     const TID& first_id,
                     const TID& last_id
                     ) :
        Author(author),
        Comment(comment),
        TimeStamp(time),
        CalibrationID(calibrationID),
        Extendable(false),
        FirstID(first_id),
        LastID(last_id),
        Data(),
        FitParameters()
    {}

    // this constructor is used by some tests
    TCalibrationData(const std::string& calibrationID,
                     const TID& first_id,
                     const TID& last_id,
                     bool extendable = false) :
        Author(),
        Comment(),
        TimeStamp(),
        CalibrationID(calibrationID),
        Extendable(extendable),
        FirstID(first_id),
        LastID(last_id),
        Data(),
        FitParameters()
    {}

    virtual std::ostream& Print( std::ostream& s) const override {
        s << "TCalibrationData:\n";
        s << "  Generated at    " << std::asctime(std::localtime(&TimeStamp)); // no << '\n'; because asctime outputs this already
        s << "  CalibrationID:  " << CalibrationID << '\n';
        s << "  Valid for IDs:  [" << FirstID << ", " << LastID << "]";
        if(Extendable)
            s << " (extendable)";
        s  << '\n';
        return s;
    }
#endif

    TCalibrationData() :
        Author(),
        Comment(),
        TimeStamp(),
        CalibrationID(),
        Extendable(false),
        FirstID(),
        LastID(),
        Data(),
        FitParameters()
    {}
    virtual ~TCalibrationData() {}
    ClassDef(TCalibrationData, ANT_CALIBRATION_DATA_VERSION)

};

}
