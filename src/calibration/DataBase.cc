#include "DataBase.h"

#include "base/WrapTFile.h"
#include "base/interval.h"
#include "base/Logger.h"

#include "base/std_ext/system.h"
#include "base/std_ext/string.h"

#include <sstream>
#include <iomanip>

using namespace std;
using namespace ant;
using namespace ant::std_ext;
using namespace ant::calibration;

DataBase::DataBase(const string calibrationDataFolder_):
    calibrationDataFolder(calibrationDataFolder_)
{

}

std::list<string> DataBase::GetCalibrationIDs() const
{
    list<string> ids;
    for(auto folder : system::lsFiles(calibrationDataFolder)) {
        auto pos_slash = folder.rfind('/');
        if(pos_slash != string::npos)
            folder = folder.substr(pos_slash+1);
        if(folder == ".." || folder == ".")
            continue;
        ids.emplace_back(folder);
    }
    return ids;
}

size_t DataBase::GetNumberOfCalibrationData(const string& calibrationID) const
{
    auto count_rootfiles = [] (const string& folder) {
        size_t n = 0;
        for(auto rootfile : system::lsFiles(folder)) {
            if(!string_ends_with(rootfile, ".root"))
                continue;
            n++;
        }
        return n;
    };

    // count the number of .root files in MC, DataDefault and DataRanges
    size_t total = 0;
    total += count_rootfiles(calibrationDataFolder+"/"+calibrationID+"/MC");
    total += count_rootfiles(calibrationDataFolder+"/"+calibrationID+"/DataDefault");
    for(auto range : getDataRanges(calibrationID)) {
        total += count_rootfiles(range.FolderPath);
    }
    return total;
}

std::set<DataBase::Range_t> DataBase::getDataRanges(const string& calibrationID) const
{
    set<Range_t> ranges;
    for(auto rangedir : system::lsFiles(calibrationDataFolder+"/"+calibrationID+"/DataRanges")) {
        /// \todo parse folder names
        cout << rangedir << endl;
    }
    return ranges;
}

bool DataBase::loadFile(const string& filename, TCalibrationData& cdata) const
{
    WrapTFileInput dataFile;
    try {
        dataFile.OpenFile(filename);
    }
    catch(...) {
        VLOG(5) << "Cannot open file for reading: " << filename;
        return false;
    }
    return dataFile.GetObjectClone("cdata", cdata);
}

bool DataBase::writeFile(const string& folder, const TCalibrationData& cdata) const
{
    // ensure the folder is there
    system::exec(formatter() << "mkdir -p " << folder);

    // get the highest numbered root file in it
    size_t maxnum = 0;
    for(auto rootfile : system::lsFiles(folder)) {
        if(!string_ends_with(rootfile, ".root"))
            continue;
        auto pos_slash = rootfile.rfind('/');
        if(pos_slash == string::npos)
            continue;
        auto pos_dot = rootfile.rfind('.');
        if(pos_dot == string::npos)
            continue;
        stringstream numstr(rootfile.substr(pos_slash+1, pos_dot-pos_slash-1));
        size_t num;
        if(numstr >> num && num>=maxnum)
            maxnum = num+1;
    }

    string filename(formatter() << setw(4) << setfill('0') << maxnum << ".root");

    try {
        WrapTFileOutput outfile(folder+"/"+filename);
        auto nbytes = outfile.WriteObject(addressof(cdata), "cdata");
        // update the symlink
        system::exec(formatter() << "cd " << folder << " && ln -s -T -f " << filename << " current");
        VLOG(5) << "Wrote TCalibrationData with " << nbytes << " bytes to " << filename;
    }
    catch(...) {
        return false;
    }

    return true;

}

bool DataBase::GetItem(const string& calibrationID,
                       const TID& currentPoint,
                       TCalibrationData& theData,
                       TID& nextChangePoint) const
{
    // handle MC
    if(currentPoint.isSet(TID::Flags_t::MC)) {
        return loadFile(calibrationDataFolder+"/"+calibrationID+"/MC/current", theData);
    }

    // try to find it in the DataRanges
    const auto ranges = getDataRanges(calibrationID);
    const auto it = ranges.find(currentPoint);
    /// \todo correctly consider ranges, watch out for default data ranges...

    // not found in ranges, so try default data
    if(loadFile(calibrationDataFolder+"/"+calibrationID+"/DataDefault/current", theData)) {
        /// \todo figure out the nextChangePoint
        VLOG(5) << "Loaded default data for " << calibrationID << " for changepoint " << currentPoint;
        return true;
    }

    // nothing found at all
    return false;
}

void DataBase::AddItem(const TCalibrationData& cdata, mode_t mode)
{
    // some checks
    if(cdata.FirstID.isSet(TID::Flags_t::MC) ^ cdata.LastID.isSet(TID::Flags_t::MC))
        throw Exception("Inconsistent flags for FirstID/LastID");

    auto& calibrationID = cdata.CalibrationID;
    if(calibrationID.empty())
        throw Exception("Provided CalibrationID is empty string");

    // handle MC
    if(cdata.FirstID.isSet(TID::Flags_t::MC))
    {
        writeFile(calibrationDataFolder+"/"+calibrationID+"/MC", cdata);
        return;
    }

    // non-MC business

    switch(mode) {
    case mode_t::AsDefault: {
        writeFile(calibrationDataFolder+"/"+calibrationID+"/DataDefault", cdata);
        break;
    }
    case mode_t::StrictRange: {
        if(cdata.FirstID.IsInvalid())
            throw Exception("FirstID cannot be invalid");
        /// \todo Implement
        break;
    }
    case mode_t::RightOpen: {
        if(cdata.FirstID.IsInvalid())
            throw Exception("FirstID cannot be invalid");
        /// \todo Implement
        break;
    }
    } // end switch

}



