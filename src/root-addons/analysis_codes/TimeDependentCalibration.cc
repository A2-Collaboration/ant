#include "TimeDependentCalibration.h"

#include "base/WrapTFile.h"
#include "analysis/plot/HistogramFactory.h"

using namespace std;
using namespace ant;

void TimeDependentCalibration::MakeCBEnergyFile(const string& outfilename)
{
    WrapTFileOutput(outfilename, WrapTFileOutput::mode_t::recreate, true);

}
