#include "CalibrationEditor.h"
#include "base/std_ext.h"

#include "TH2D.h"


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


void Editor::ShowHistory(const string& calibrationID) const
{
    auto& dVector = dman.DataMap.at(calibrationID);

    interval<TID> maxInt(TID(0,0),TID(0,0));
    if (!getIDRange(calibrationID,maxInt))
    {
        cout << "Calibration Doesn't exist" << endl;
        return;
    }
    auto len = maxInt.Stop().Value - maxInt.Start().Value;
    auto first = maxInt.Start().Value;

    TH2D* hist = new TH2D( (std_ext::formatter() << "hist-" << calibrationID).str().c_str(),
                          (std_ext::formatter() << "History for " << calibrationID).str().c_str(),
                          100, 0, 100,
                          dman.GetNumberOfDataPoints(calibrationID),
                          0,
                          dman.GetNumberOfDataPoints(calibrationID)
                          );

    uint32_t step(0);
    for(const auto& cdata: dVector)
    {
        unsigned min = (unsigned) (cdata.FirstID.Value - first) * 100 / len;
        unsigned max = (unsigned) (cdata.LastID.Value - first) * 100 / len;

        for ( auto i = min ; i <= max ; ++i)
        {
            hist->Fill(i,step);
        }
        ++step;
        cout << "show result" << endl;
    }

}
