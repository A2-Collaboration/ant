#pragma once

namespace ant {
struct TimeDependentCalibration {
    static void MakeCBEnergyFile(const char* basefilename,
                                 const char* setupname,
                                 int fillsPerChannel = 10000,
                                 int nSlices = 10);
};

}
