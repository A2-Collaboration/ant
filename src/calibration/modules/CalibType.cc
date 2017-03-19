#include "CalibType.h"

using namespace std;
using namespace ant;
using namespace ant::calibration;

double CalibType::Get(unsigned channel) const {

    if(Values.empty()) {
        if(DefaultValues.size() == 1) {
            return DefaultValues.front();
        }
        else {
            return DefaultValues.at(channel);
        }
    }
    else {
        return Values.at(channel);
    }
}

CalibType::CalibType(
        const std::shared_ptr<const Detector_t>& det,
        const string& name,
        const std::vector<double>& defaultValues,
        const string& histname) :
    Name(name),
    // use name for histogram if not provided different
    HistogramName(histname.empty() ? name : histname),
    Values(),
    DefaultValues(defaultValues)
{
    if(DefaultValues.size() != 1 && DefaultValues.size() != det->GetNChannels()) {
        throw runtime_error("Wrong size of default values for calibType="+name+" det="+Detector_t::ToString(det->Type));
    }
}