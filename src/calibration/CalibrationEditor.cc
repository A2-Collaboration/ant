#include "CalibrationEditor.h"

using namespace std;
using namespace ant::calibration;

bool Editor::getIDRange(const string& calibrationID, interval<TID>& IDinterval) const
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

void Editor::PrintHistory(const string& calibrationID) const
{
    auto& dVector = dman.DataMap.at(calibrationID);

    interval<TID> maxInt(TID(0,0),TID(0,0));
    if (!getIDRange(calibrationID,maxInt))
    {
        cout << "Calibration Doesn't exist" << endl;
        return;
    }



    for(auto& cdata: dVector)
    {
        cout << "show result" << endl;
    }

}
