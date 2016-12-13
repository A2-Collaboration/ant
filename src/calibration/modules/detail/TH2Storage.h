#pragma once

#include "base/BinSettings.h"

class TH2;
class TH2D;

namespace ant {

class TCalibrationData;

namespace calibration {
namespace detail {

struct TH2Storage {
    static void Encode(const TH2* histogram, TCalibrationData& cdata);
    static TH2D* Decode(const TCalibrationData& cdata);
    static BinSettings GetXBins(const TCalibrationData& cdata);
    static BinSettings GetYBins(const TCalibrationData& cdata);
};

}
}
}
