#include "DataBase.h"

#include "base/WrapTFile.h"
#include "base/interval.h"
#include "base/Logger.h"

#include "base/std_ext/system.h"
#include "base/std_ext/string.h"
#include "base/std_ext/time.h"
#include "base/std_ext/misc.h"
#include "base/std_ext/math.h"


#include <sstream>
#include <iomanip>
#include <ctime>

using namespace std;
using namespace ant;
using namespace ant::std_ext;
using namespace ant::calibration;

DataBase::DataBase(const string& calibrationDataFolder):
    Layout(calibrationDataFolder)
{

}

bool DataBase::GetItem(const string& calibrationID,
                       const TID& currentPoint,
                       TCalibrationData& theData,
                       TID& nextChangePoint) const
{
    nextChangePoint = TID();

    // handle MC
    if(currentPoint.isSet(TID::Flags_t::MC)) {
        if(loadFile(Layout.GetCurrentFile(calibrationID, OnDiskLayout::Type_t::MC), theData)) {
            LOG(INFO) << "Loaded MC data for " << calibrationID;
            return true;
        }
        return false;
    }

    // try to find it in the DataRanges
    // the following needs the ranges to be sorted after their
    auto ranges = Layout.GetDataRanges(calibrationID);

    auto it_range = find_if(ranges.begin(), ranges.end(),
                            [currentPoint] (const OnDiskLayout::Range_t& r) {
        if(r.Stop().IsInvalid())
            return r.Start() < currentPoint;
        return r.Contains(currentPoint);
    });

    if(it_range != ranges.end()) {
        if(loadFile(it_range->FolderPath+"/current", theData)) {
            LOG(INFO) << "Loaded data for " << calibrationID << " for changepoint " << currentPoint << " from " << it_range->FolderPath;
            // next change point is given by found range
            nextChangePoint = it_range->Stop();
            ++nextChangePoint;
            return true;
        }
        else {
            LOG(WARNING) << "Cannot load data from " << it_range->FolderPath;
        }
    }

    // check if there's a range coming up at some point
    // that means even if this method returns false,
    // the nextChangePoint is correctly set
    ranges.sort();
    it_range = find_if(ranges.begin(), ranges.end(),
                       [currentPoint] (const OnDiskLayout::Range_t& r) {
        return currentPoint < r.Start();
    });
    if(it_range != ranges.end())
        nextChangePoint = it_range->Start();

    // not found in ranges, so try default data
    if(loadFile(Layout.GetCurrentFile(calibrationID, OnDiskLayout::Type_t::DataDefault), theData)) {
        LOG(INFO) << "Loaded default data for " << calibrationID << " for changepoint " << currentPoint;
        return true;
    }

    // nothing found at all
    return false;
}

void DataBase::AddItem(const TCalibrationData& cdata, Calibration::AddMode_t mode)
{
    // some general checks
    if(cdata.FirstID.isSet(TID::Flags_t::MC) ^ cdata.LastID.isSet(TID::Flags_t::MC))
        throw Exception("Inconsistent flags for FirstID/LastID");

    auto& calibrationID = cdata.CalibrationID;
    if(calibrationID.empty())
        throw Exception("Provided CalibrationID is empty string");

    // handle MC
    if(cdata.FirstID.isSet(TID::Flags_t::MC))
    {
        writeToFolder(Layout.GetFolder(calibrationID, OnDiskLayout::Type_t::MC), cdata);
        return;
    }

    // non-MC business

    switch(mode) {
    case Calibration::AddMode_t::AsDefault: {
        writeToFolder(Layout.GetFolder(calibrationID, OnDiskLayout::Type_t::DataDefault), cdata);
        break;
    }
    case Calibration::AddMode_t::StrictRange: {
        handleStrictRange(cdata);
        break;
    }
    case Calibration::AddMode_t::RightOpen: {
        handleRightOpen(cdata);
        break;
    }
    } // end switch

}

void DataBase::handleStrictRange(const TCalibrationData& cdata) const
{
    if(cdata.FirstID.IsInvalid())
        throw Exception("FirstID cannot be invalid");
    if(cdata.LastID.IsInvalid())
        throw Exception("LastID cannot be invalid");
    if(cdata.LastID < cdata.FirstID)
        throw Exception("LastID<FirstID cannot be");

    const auto& calibrationID = cdata.CalibrationID;
    const interval<TID> range(cdata.FirstID, cdata.LastID);

    // check if range already exists (ranges don't need to be sorted for this)
    const auto ranges = Layout.GetDataRanges(calibrationID);
    const auto it_range = find_if(ranges.begin(), ranges.end(),
                                  [range] (const OnDiskLayout::Range_t& r) {
        return !r.Disjoint(range);
    });

    // given range is disjoint with all existing ones,
    // then just add this range
    if(it_range != ranges.end()) {
        // if range overlaps, then they must exactly match
        // anything else is not allowed right now
        if(*it_range != range) {
            if(std_ext::abs_diff(it_range->Stop().Lower, range.Stop().Lower) < 100) {
                LOG(WARNING) << "Fixing TCalibrationData range: " << *it_range  << " to " << range;

                const auto oldfolder = Layout.GetRangeFolder(calibrationID, *it_range);
                const auto newfolder = Layout.GetRangeFolder(calibrationID, range);

                system::exec(formatter() << "mv " << oldfolder << " " << newfolder);

            } else
                throw Exception(formatter() << "Given TCalibrationData range(" << range << ") conflicts with existing database entry(" << *it_range << ").");
        }
    }

    writeToFolder(Layout.GetRangeFolder(calibrationID, range), cdata);
}

void DataBase::handleRightOpen(const TCalibrationData& cdata) const
{
    if(cdata.FirstID.IsInvalid())
        throw Exception("FirstID cannot be invalid");

    const auto& calibrationID = cdata.CalibrationID;
    // LastID of cdata is simply ignored
    const auto startPoint = cdata.FirstID;
    interval<TID> range(startPoint, TID());

    // scan the ranges for conflicts
    auto ranges = Layout.GetDataRanges(calibrationID);
    ranges.sort();
    auto it_conflict = find_if(ranges.begin(), ranges.end(),
                            [startPoint] (const OnDiskLayout::Range_t& r) {
        // two half-open intervals always conflict
        if(r.Stop().IsInvalid())
            return true;
        // conflict if startPoint is before or equal the stop of existing interval
        return startPoint <= r.Stop();
    });

    if(it_conflict != ranges.end()) {

        // handle conflicts depending on where the given startPoint is

        if(startPoint == it_conflict->Start()) {
            // add to existing
           range = *it_conflict;
        }
        else if(startPoint > it_conflict->Start()) {
            // shrink existing by renaming folder
            it_conflict->Stop() = startPoint;
            --(it_conflict->Stop());
            const auto& newfolder = Layout.GetRangeFolder(calibrationID, *it_conflict);
            system::exec(formatter() << "mv " << it_conflict->FolderPath
                         << " " << newfolder);

            // search if new data has some stop
            // important that ranges were sorted by start
            ++it_conflict;
            if(it_conflict != ranges.end()) {
                range.Stop() = it_conflict->Start();
                --range.Stop();
            }
        }
        else if(startPoint < it_conflict->Start()) {
            // shrink new data
            range.Stop() = it_conflict->Start();
            --range.Stop();
        }

    }

    writeToFolder(Layout.GetRangeFolder(calibrationID, range), cdata);
}

std::list<string> DataBase::GetCalibrationIDs() const
{
    return system::lsFiles(Layout.CalibrationDataFolder,"",true,true);
}

size_t DataBase::GetNumberOfCalibrationData(const string& calibrationID) const
{
    auto count_rootfiles = [] (const string& folder) {
        return system::lsFiles(folder, ".root").size();
    };

    // count the number of .root files in MC, DataDefault and DataRanges
    size_t total = 0;
    total += count_rootfiles(Layout.GetFolder(calibrationID, OnDiskLayout::Type_t::MC));
    total += count_rootfiles(Layout.GetFolder(calibrationID, OnDiskLayout::Type_t::DataDefault));
    for(auto range : Layout.GetDataRanges(calibrationID)) {
        total += count_rootfiles(range.FolderPath);
    }
    return total;
}

bool DataBase::loadFile(const string& filename, TCalibrationData& cdata) const
{
    string errmsg;
    if(!system::testopen(filename, errmsg)) {
        VLOG(8) << "Cannot open " << filename << ": " << errmsg;
        return false;
    }


    try {
        WrapTFileInput dataFile;
        dataFile.OpenFile(filename);
        return dataFile.GetObjectClone("cdata", cdata);
    }
    catch(...) {
        LOG(WARNING) << "Cannot load object cdata from " << filename;
        return false;
    }
}

bool DataBase::writeToFolder(const string& folder, const TCalibrationData& cdata) const
{
    // ensure the folder is there
    system::exec(formatter() << "mkdir -p " << folder);

    // get the highest numbered root file in it
    size_t maxnum = 0;
    for(auto rootfile : system::lsFiles(folder, ".root",true,true)) {
        auto pos_dot = rootfile.rfind('.');
        if(pos_dot == string::npos)
            continue;
        stringstream numstr(rootfile.substr(0, pos_dot));
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

string DataBase::OnDiskLayout::GetRangeFolder(const string& calibrationID, const interval<TID>& range) const
{
    const string& start = makeTIDString(range.Start());
    const string& stop = makeTIDString(range.Stop());
    const string& day = start.substr(0, start.find('T'));
    return GetFolder(calibrationID, Type_t::DataRanges)+"/"+day+"/"+start+"-"+stop;
}

string DataBase::OnDiskLayout::GetCurrentFile(const DataBase::OnDiskLayout::Range_t& range) const
{
    return range.FolderPath + "/current";
}

std::list<DataBase::OnDiskLayout::Range_t> DataBase::OnDiskLayout::GetDataRanges(const string& calibrationID) const
{

    auto it_cached_range = cached_ranges.find(calibrationID);
    if(it_cached_range != cached_ranges.end())
        return it_cached_range->second;

    list<Range_t> ranges;
    for(auto daydir : system::lsFiles(GetFolder(calibrationID, Type_t::DataRanges),"",true)) {
        for(auto tidRangeDir : system::lsFiles(daydir, "", true, true)) {
            auto tidRange = parseTIDRange(tidRangeDir);
            ranges.emplace_back(tidRange, daydir+"/"+tidRangeDir);
        }
    }

    if(EnableCaching)
        cached_ranges.emplace(calibrationID, ranges);

    return ranges;
}

bool DataBase::OnDiskLayout::EnableCaching = false;

DataBase::OnDiskLayout::OnDiskLayout(const string& calibrationDataFolder) :
    CalibrationDataFolder(calibrationDataFolder) {}

string DataBase::OnDiskLayout::GetFolder(const string& calibrationID, DataBase::OnDiskLayout::Type_t type) const
{
    switch(type) {
    case Type_t::DataDefault: return CalibrationDataFolder+"/"+calibrationID+"/DataDefault";
    case Type_t::DataRanges: return CalibrationDataFolder+"/"+calibrationID+"/DataRanges";
    case Type_t::MC: return CalibrationDataFolder+"/"+calibrationID+"/MC";
    }
    throw runtime_error("Invalid folder type provided.");
}

string DataBase::OnDiskLayout::GetCurrentFile(const string& calibrationID, DataBase::OnDiskLayout::Type_t type) const
{
    switch(type) {
    case Type_t::DataDefault:
    case Type_t::MC:
        return GetFolder(calibrationID, type)+"/current";
    default:
        throw runtime_error("Invalid folder type provided.");
    }

}

string DataBase::OnDiskLayout::makeTIDString(const TID& tid) const
{
    if(tid.IsInvalid())
        return "OPEN";
    // cannot use std_ext::to_iso8601 since
    // only underscores are allowed in string
    // note that the T is important for some use cases
    char buf[sizeof "2011_10_08T07_07_09Z"];
    time_t time = tid.Timestamp;
    strftime(buf, sizeof buf, "%Y_%m_%dT%H_%M_%SZ", gmtime(addressof(time)));
    // add the lower ID
    stringstream ss;
    ss << buf << "_0x" << hex << setw(8) << setfill('0') << tid.Lower;
    return ss.str();
}

interval<TID> DataBase::OnDiskLayout::parseTIDRange(const string& tidRangeStr) const
{
    auto pos_dash = tidRangeStr.find('-');
    if(pos_dash == string::npos)
        return interval<TID>(TID(), TID());
    auto startStr = tidRangeStr.substr(0,pos_dash);
    auto stopStr = tidRangeStr.substr(pos_dash+1);

    auto parseTIDStr = [] (const string& str) {
        if(str == "OPEN")
            return TID();
        if(str.size() != 31)
            return TID();
        if(str.substr(20,3) != "_0x")
            return TID();

        // the following calculation assumes
        // UTC as timezone
        const char* tz = getenv("TZ");
        setenv("TZ", "UTC", 1);
        tzset();

        std_ext::execute_on_destroy restore_TZ([tz] () {
            // restore TZ environment
            if(tz)
                setenv("TZ", tz, 1);
            else
                unsetenv("TZ");
            tzset();
        });
        auto timestr = str.substr(0, 20);
        time_t timestamp = to_time_t(to_tm(timestr, "%Y_%m_%dT%H_%M_%SZ"));

        stringstream lowerstr;
        lowerstr << hex << str.substr(23); // skip leading 0x
        std::uint32_t lower;
        if(!(lowerstr >> lower))
            return TID();

        return TID(timestamp, lower);
    };

    return {parseTIDStr(startStr), parseTIDStr(stopStr)};
}
