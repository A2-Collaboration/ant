#include "DataManager.h"

//ant
#include "base/WrapTFile.h"
#include "base/Logger.h"
#include "base/std_ext/memory.h"
#include "base/interval.h"
#include "DataBase.h"
#include "tree/TCalibrationData.h"
#include "tree/TDataRecord.h"

//ROOT
#include "TTree.h"
#include "TError.h"
#include "TDirectory.h"

//std
#include <algorithm>
#include <iostream>

using namespace std;
using namespace ant;
using namespace ant::calibration;


void DataManager::Init()
{
    if(dataBase)
        return;
    dataBase = std_ext::make_unique<DataBase>();
    dataBase->ReadFromFolder(calibrationDataFolder);
}

DataManager::DataManager(const string& calibrationDataFolder_):
    calibrationDataFolder(calibrationDataFolder_),
    extendable(false)
{}

DataManager::~DataManager()
{
    if(dataBase)
        dataBase->WriteToFolder(calibrationDataFolder);
}

void DataManager::Add(const TCalibrationData& cdata)
{
    Init();
    // explicitly copy the given object
    // to modify its extendable flag
    auto cdata_ = cdata;
    if(extendable) {
        VLOG(3) << "Flagged given TCalibrationData as extendable";
        cdata_.Extendable = true;
    }
    dataBase->AddItem(cdata_);
    LOG(INFO) << "Added " << cdata_;
}


bool DataManager::GetData(const string& calibrationID,
                          const TID& eventID, TCalibrationData& cdata)
{
    Init();
    // case one: calibration doesn't exist at all
    if (!dataBase->Has(calibrationID))
        return false;

    // case two: eventID lies inside data, then return it
    // but also track TCalibrationData items which are extendable

    using extendable_item_t = pair<TID, const TCalibrationData*>;
    vector<extendable_item_t> forwards;
    vector<extendable_item_t> backwards;

    const auto& calibPairs = dataBase->GetItems(calibrationID);
    for(auto rit = calibPairs.rbegin(); rit != calibPairs.rend(); ++rit)
    {
        interval<TID> range(rit->FirstID, rit->LastID);
        if (range.Contains(eventID))
        {
            cdata = *rit;
            return true;
        }
        const bool extendable = rit->Extendable
                                || (rit->FirstID.isSet(TID::Flags_t::MC) && rit->FirstID.isSet(TID::Flags_t::MC));
        if(extendable
           && rit->FirstID.Flags == eventID.Flags
           && rit->LastID.Flags == eventID.Flags)
        {
            forwards.emplace_back(rit->FirstID, addressof(*rit));
            backwards.emplace_back(rit->LastID, addressof(*rit));
        }
    }

    // case three: there were extendable calibration data
    if(!forwards.empty() && !backwards.empty()) {
        // use stable_sort in case earlier added TCalibrationData has same TID start/stop
        auto forwards_comparer = [] (const extendable_item_t& a, const extendable_item_t& b) {
            return a.first < b.first;
        };
        auto backwards_comparer = [] (const extendable_item_t& a, const extendable_item_t& b) {
            return a.first > b.first;
        };
        stable_sort(forwards.begin(), forwards.end(), forwards_comparer);
        stable_sort(backwards.begin(), backwards.end(), backwards_comparer);

        // search position of given eventID
        auto it_forward = forwards.begin();
        while(it_forward != forwards.end() && eventID > it_forward->first)
            ++it_forward;
        if(it_forward != forwards.end()) {
            if(eventID > it_forward->first)
                --it_forward;
            cdata = *(it_forward->second);
            return true;
        }

        auto it_backward = backwards.begin();
        while(it_backward != backwards.end() && eventID < it_backward->first)
            ++it_backward;
        if(it_backward != backwards.end()) {
            if(eventID < it_backward->first)
                --it_backward;
            cdata = *(it_backward->second);
            return true;
        }

        throw runtime_error("Should be never reached");
    }


    // case four: TID not covered by anything
    return false;
}


const list<TID> DataManager::GetChangePoints(const string& calibrationID)
{
    Init();
    return dataBase->GetChangePoints(calibrationID);
}

uint32_t DataManager::GetNumberOfCalibrations()
{
    Init();
    return dataBase->GetKeys().size();
}

uint32_t DataManager::GetNumberOfDataPoints(const string& calibrationID)
{
    Init();
    return dataBase->GetNumberOfDataPoints(calibrationID);
}
