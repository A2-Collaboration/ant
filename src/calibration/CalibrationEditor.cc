#include "CalibrationEditor.h"

using namespace std;
using namespace ant::calibration;

bool Editor::getIDRange(const string& calibrationID, interval<TID>& IDinterval)
{
    if (dman.DataMap.count(calibrationID) == 0)
        return false;

    auto& data = dman.DataMap.at(calibrationID);

    IDinterval.Start() = data.front().FirstID;
    IDinterval.Stop()  = data.front().LastID;

    for (auto& entry: data)
    {
        if (entry.FirstID < IDinterval.Start())
                IDinterval.Start() = entry.FirstID;
        if (entry.LastID  > IDinterval.Stop())
                IDinterval.Stop() = entry.LastID;
    }

    return true;
}
