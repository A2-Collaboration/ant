#pragma once

#include "TID.h"
#include "base/GitInfo.h"

#include <string>
#include <ctime>

#define ANT_CALIBRATION_DATA_VERSION 1

namespace ant {

/**
 * @brief The TCalibrationData struct contains one calibration iteration
 */
struct TCalibrationData
{
    std::string Author;

    std::int64_t TimeStamp;

    std::string CalibrationID;

    bool Extendable; // OBSOLETE
    TID FirstID;
    TID LastID;

    typedef TKeyValue<double> Entry;
    std::vector<Entry> Data;

    typedef TKeyValue<std::vector<double>> TFitParameters;
    std::vector<TFitParameters> FitParameters;

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

    friend std::ostream& operator<<( std::ostream& s, const TCalibrationData& o) {
        s << "TCalibrationData:\n";
        s << "  Generated at    " << std_ext::to_iso8601(o.TimeStamp) << " by " << o.Author << '\n';
        s << "  CalibrationID:  " << o.CalibrationID << '\n';
        s << "  Valid for IDs:  [" << o.FirstID << ", " << o.LastID << "]\n";

        double avgFitParameters = 0;
        for(const TFitParameters& fitparam : o.FitParameters)
            avgFitParameters += fitparam.Value.size();
        avgFitParameters /= o.FitParameters.size();
        s << "  nEntries=" << o.Data.size() << " nFitParameters=" << o.FitParameters.size()
          << " avgFitParameters=" << avgFitParameters << '\n';

        return s;
    }

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
