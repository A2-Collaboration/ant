#pragma once

#include <string>

class TH1;
class TH2;
class TF1;

namespace ant {
class TimeDependentCalibration {
    static void MakeCBEnergyFile(const std::string& outfilename);
};

}
