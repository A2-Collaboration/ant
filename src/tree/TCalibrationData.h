#pragma once

#include "TID.h"

#include <string>
#include <ctime>

#ifndef __CINT__
#include "base/GitInfo.h"
#endif

#define ANT_CALIBRATION_DATA_VERSION 1

namespace ant {

/**
 * @brief The TCalibrationData struct contains one calibration iteration
 */
#ifndef __CINT__
struct TCalibrationData : printable_traits
#else
struct TCalibrationData
#endif
{
    std::string Author;

    std::int64_t TimeStamp;

    std::string CalibrationID;

    bool Extendable; // OBSOLETE
    TID FirstID;
    TID LastID;

    typedef TKeyValue<double> TCalibrationEntry;
    std::vector<TCalibrationEntry> Data;

    typedef TKeyValue<std::vector<double>> TFitParameters;
    std::vector<TFitParameters> FitParameters;

#ifndef __CINT__
    TCalibrationData(const std::string& calibrationID,
                     const TID& first_id,
                     const TID& last_id
                     ) :
        Author(),
        TimeStamp(std::time(nullptr)),
        CalibrationID(calibrationID),
        FirstID(first_id),
        LastID(last_id),
        Data(),
        FitParameters()
    {
        GitInfo info;
        Author = info.GetUser();
    }

    virtual std::ostream& Print( std::ostream& s) const override {
        s << "TCalibrationData:\n";
        s << "  Generated at    " << std_ext::to_iso8601(TimeStamp) << " by " << Author << '\n';
        s << "  CalibrationID:  " << CalibrationID << '\n';
        s << "  Valid for IDs:  [" << FirstID << ", " << LastID << "]";
        s  << '\n';
        return s;
    }
#endif

    TCalibrationData() :
        Author(),
        TimeStamp(),
        CalibrationID(),
        FirstID(),
        LastID(),
        Data(),
        FitParameters()
    {}
    virtual ~TCalibrationData() {}
    ClassDef(TCalibrationData, ANT_CALIBRATION_DATA_VERSION)

};

}
