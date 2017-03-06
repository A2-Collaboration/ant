#pragma once

namespace ant {
struct TimeDependentCalibration {
    static void MakeCBEnergyFile(const char* outfilename,
                                 int fillsPerChannel = 10000);
};

}
