#include "CalibrationEditor.h"
#include "base/std_ext.h"
#include "tree/TCalibrationData.h"

#include "TH2D.h"
#include "algorithm"


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

vector<pair<uint32_t,TCalibrationData>> Editor::getValidData(const std::string& calibrationID) const
{
    vector<pair<uint32_t,TCalibrationData>> theList;
    if (dman.DataMap.count(calibrationID) == 0)
        return theList;
    cout << "HERE I AM:" << endl;

    auto changePoints = dman.GetChangePoints(calibrationID);
    auto& dVector = dman.DataMap.at(calibrationID);
    cout << dVector.size() << endl;
    cout << changePoints.size() << endl;

    for (const auto&  cp: changePoints)
    {
        cout << cp.Value << endl;
        uint32_t index = dVector.size();
        for(auto rit = dVector.rbegin(); rit != dVector.rend(); ++rit)
        {
            --index;
            interval<TID> rit_int(rit->FirstID,rit->LastID);
            if (rit_int.Contains(cp))
            {
                pair<uint32_t,TCalibrationData> thePair;
                thePair.first = index;
                thePair.second = *rit;
                theList.push_back(thePair);
                break;
            }
        }
    }


    //remove doubles
    std::sort(theList.begin(),theList.end(),[](const pair<uint32_t,TCalibrationData>& a,
                                               const pair<uint32_t,TCalibrationData>& b) -> bool
                                               {
                                                   return a.first < b.first;
                                               } );
    auto last = std::unique(theList.begin(),theList.end(),[](const pair<uint32_t,TCalibrationData>& a,
                                               const pair<uint32_t,TCalibrationData>& b) -> bool
                                               {
                                                   return a.first == b.first;
                                               });
    theList.erase(last,theList.end());

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

    TH2D* valid = new TH2D( (std_ext::formatter() << "val-" << calibrationID).str().c_str(),
                          (std_ext::formatter() << "Valid data for " << calibrationID).str().c_str(),
                          100, 0, 100,
                            dman.GetNumberOfDataPoints(calibrationID),
                          0,
                            dman.GetNumberOfDataPoints(calibrationID)
                          );
    valid->SetXTitle("TID [%]");
    valid->SetYTitle("crap, has no meaning");

    for(const auto& cdata: dVector)
    {
        unsigned min = (unsigned) (cdata.second.FirstID.Value - first) * 100 / len;
        unsigned max = (unsigned) (cdata.second.LastID.Value - first) * 100 / len;

        for ( auto i = min ; i <= max ; ++i)
            valid->Fill(i,cdata.first);
    }

    valid->Draw("col");

}
