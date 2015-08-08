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


bool Editor::ShowHistory(const string& calibrationID) const
{

    interval<TID> maxInt(TID(0,0),TID(0,0));
    if (!getIDRange(calibrationID,maxInt))
        return false;

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


    //Fill ranges of calibration data
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
    //Fill valid IDs twice
    auto valids = getValidData(calibrationID);
    for (const auto& vcdata: valids)
    {
        unsigned min = (unsigned) (vcdata.second.FirstID.Value - first) * 100 / len;
        unsigned max = (unsigned) (vcdata.second.LastID.Value - first) * 100 / len;

        for ( auto i = min ; i <= max ; ++i)
        {
            hist->Fill(i,vcdata.first);
        }

    }

    hist->Draw("col");

    return true;

}

vector<pair<uint32_t,TCalibrationData>> Editor::getValidData(const std::string& calibrationID) const
{
    vector<pair<uint32_t,TCalibrationData>> theList;

    if ( !dman.Has(calibrationID) )
        return theList;

    auto changePoints = dman.GetChangePoints(calibrationID);
    auto& dVector = dman.DataMap.at(calibrationID);

    for (const auto&  cp: changePoints)
    {
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


    //sort
    std::sort(theList.begin(),theList.end(),
              [](const pair<uint32_t,TCalibrationData>& a,
                 const pair<uint32_t,TCalibrationData>& b) -> bool
                 {
                     return a.first < b.first;
                 } );
    //remove doubles
    auto last = std::unique(theList.begin(),theList.end(),
                            [](const pair<uint32_t,TCalibrationData>& a,
                               const pair<uint32_t,TCalibrationData>& b) -> bool
                               {
                                   return a.first == b.first;
                               } );
    theList.erase(last,theList.end());

    return theList;
}

bool Editor::Remove(const string &calibrationID, const uint32_t &index)
{
    if (!dman.Has(calibrationID))
        return false;

    auto& dVector = dman.DataMap.at(calibrationID);
    if (index >= dVector.size())
        return false;
    dVector.erase(dVector.begin()+index);
    return true;
}

bool Editor::Remove(const string &calibrationID, const uint32_t &index1, const uint32_t &index2)
{
    if (!dman.Has(calibrationID))
        return false;

    auto& dVector = dman.DataMap.at(calibrationID);
    if ((index1 >= dVector.size()) && (index2 >= dVector.size()))
        return false;
    if ( index1 > index2 )
        return false;

    dVector.erase(dVector.begin()+index1,dVector.begin()+index2);
    return true;

}


bool Editor::ShowValid(const string& calibrationID) const
{
    interval<TID> maxInt(TID(0,0),TID(0,0));
    if (!getIDRange(calibrationID,maxInt))
    {
        return false;
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
    valid->SetYTitle("step");

    for(const auto& cdata: dVector)
    {
        unsigned min = (unsigned) (cdata.second.FirstID.Value - first) * 100 / len;
        unsigned max = (unsigned) (cdata.second.LastID.Value - first) * 100 / len;

        for ( auto i = min ; i <= max ; ++i)
            valid->Fill(i,cdata.first);
    }

    valid->Draw("col");

    return true;
}

bool Editor::ReduceToValid(const string &calibrationID)
{
   if(!dman.Has(calibrationID))
       return false;

   auto valids = getValidData(calibrationID);

   dman.DataMap.at(calibrationID).clear();
   for (const auto& pair: valids)
       Add(pair.second);

   return true;
}
