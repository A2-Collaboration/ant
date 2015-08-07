#include "CalibrationEditor.h"
#include "base/std_ext.h"
#include "tree/TCalibrationData.h"

#include "TH2D.h"


using namespace std;
using namespace ant;
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

void Editor::ListCalibrations() const
{
    cout << "    Available Calibrations:" << endl << endl;
    for (const auto& entry: dman.DataMap){
        cout << entry.first << endl;
    }
    cout << endl;
}


void Editor::ShowHistory(const string& calibrationID) const
{

    interval<TID> maxInt(TID(0,0),TID(0,0));
    if (!getIDRange(calibrationID,maxInt))
    {
        cout << "Calibration Doesn't exist" << endl;
        return;
    }
    auto& dVector = dman.DataMap.at(calibrationID);
    auto len = maxInt.Stop().Value - maxInt.Start().Value;
    auto first = maxInt.Start().Value;

    TH2D* hist = new TH2D( (std_ext::formatter() << "hist-" << calibrationID).str().c_str(),
                          (std_ext::formatter() << "History for " << calibrationID).str().c_str(),
                          100, 0, 100,
                          dman.GetNumberOfDataPoints(calibrationID),
                          0,
                          dman.GetNumberOfDataPoints(calibrationID)
                          );
    hist->SetXTitle("TID [%]");
    hist->SetYTitle("Calibration Step");

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
    }
    hist->Draw("col");

}

vector<TCalibrationData> Editor::getValidData(const std::string& calibrationID) const
{
    vector<TCalibrationData> theList;
    interval<TID> maxInt(TID(0,0),TID(0,0));
    if (!getIDRange(calibrationID,maxInt))
        return theList;
    for(auto rit = dman.DataMap.at(calibrationID).rbegin();
        rit != dman.DataMap.at(calibrationID).rend();
        ++rit)
    {
        interval<TID> currentInt(rit->FirstID,rit->LastID);
        bool accepted = true;
        for (auto& accCal: theList){
            interval<TID> accInt(accCal.FirstID,accCal.LastID);
            if (!currentInt.Disjoint(accInt))
            {
                accepted = false;
                break;
            }
        }
        if (accepted)
            theList.push_back(*rit);

    }
    return theList;
}


void Editor::ShowValid(const string& calibrationID) const
{
    interval<TID> maxInt(TID(0,0),TID(0,0));
    if (!getIDRange(calibrationID,maxInt))
    {
        cout << "Calibration Doesn't exist" << endl;
        return;
    }
    auto dVector = getValidData(calibrationID);
    auto len = maxInt.Stop().Value - maxInt.Start().Value;
    auto first = maxInt.Start().Value;

    TH2D* hist = new TH2D( (std_ext::formatter() << "val-" << calibrationID).str().c_str(),
                          (std_ext::formatter() << "Valid data for " << calibrationID).str().c_str(),
                          100, 0, 100,
                           dVector.size(),
                          0,
                           dVector.size()
                          );
    hist->SetXTitle("TID [%]");
    hist->SetYTitle("crap, has no meaning");

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
    }
    hist->Draw("col");

}
