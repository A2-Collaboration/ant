#include "calibration/DataManager.h"
#include "calibration/DataBase.h"
#include "calibration/gui/AvgBuffer.h"

#include "expconfig/ExpConfig.h"

#include "tree/TCalibrationData.h"

#include "base/CmdLine.h"
#include "base/detail/tclap/ValuesConstraintExtra.h"
#include "base/std_ext/system.h"
#include "base/piecewise_interval.h"
#include "base/Logger.h"

#include "detail/tools.h"

using namespace std;
using namespace ant;
using namespace ant::calibration;

namespace ant {
namespace calibration {
namespace gui {
// provide AvgBuffer traits for SavitzyGolay, where Add() method is not needed
template<>
struct AvgBufferItem_traits<TCalibrationData> {
    static std::unique_ptr<TCalibrationData> Clone(const TCalibrationData& cdata) {
        // just use copy ctor
        return std_ext::make_unique<TCalibrationData>(cdata);
    }
    static int    GetNBins(const TCalibrationData& cdata) { return cdata.Data.size(); };
    static double GetBin(const TCalibrationData& cdata, int bin) { return cdata.Data[bin].Value; }
    static void   SetBin(TCalibrationData& cdata, int bin, double v) { cdata.Data[bin].Value = v; }
};
}}} // namespace ant::calibration::gui


// ensure same mapping of index in cdata.Data to Key
bool check_compatibility(const TCalibrationData& cdata) {
    static TCalibrationData prev_cdata;
    if(prev_cdata.Data.empty()) {
        prev_cdata = cdata;
        // do accept everything except empty TCalibrationData.Data
        return !cdata.Data.empty();
    }
    if(cdata.Data.size() != prev_cdata.Data.size())
        return false;
    for(auto i=0u;i<cdata.Data.size();i++) {
        if(cdata.Data[i].Key != prev_cdata.Data[i].Key)
            return false;
    }
    return true;
}

int main(int argc, char** argv)
{
    SetupLogger();

    TCLAP::CmdLine cmd("Ant-calib-smooth - applies low-pass filter on calibration data ranges", ' ', "0.1");

    TCLAP::ValuesConstraintExtra<decltype(ExpConfig::Setup::GetNames())> allowedsetupnames(ExpConfig::Setup::GetNames());
    auto cmd_setup  = cmd.add<TCLAP::ValueArg<string>>("s","setup","Use setup to determine calibration database path",true,"", &allowedsetupnames);
    auto cmd_calibration = cmd.add<TCLAP::ValueArg<string>>("c","calibration","Calibration ID", true, "","calibration");
    auto cmd_average = cmd.add<TCLAP::ValueArg<unsigned>>("a","average","Average length for Savitzky-Golay filter", false, 10, "length");
    auto cmd_sgpol = cmd.add<TCLAP::ValueArg<unsigned>>("","polyorder","Polynom order for Savitzky-Golay filter (zero is moving average)", false, 4, "polorder");
    auto cmd_channels = cmd.add<TCLAP::MultiArg<string>>("","ch","Ranges of channels, e.g. 400-412",false,"channels");
    auto cmd_dump  = cmd.add<TCLAP::SwitchArg>("","dump","Dump to stdout",false);
    auto cmd_write  = cmd.add<TCLAP::SwitchArg>("","write","Write to database",false);

    cmd.parse(argc, argv);

    if(!cmd_dump->isSet() && cmd_channels->isSet()) {
        LOG(ERROR) << "Using --channels with --dump makes no sense, as --write writes all to database";
        return EXIT_FAILURE;
    }

    if(!(cmd_dump->isSet() ^ cmd_write->isSet())) {
        LOG(ERROR) << "Use either --dump or --write, but not both.";
        return EXIT_FAILURE;
    }

    const auto channels = progs::tools::parse_cmdline_ranges(cmd_channels->getValue());
    const auto dump = cmd_dump->getValue();
    const auto write = cmd_write->getValue();

    // enable caching of the calibration database
    DataBase::OnDiskLayout::EnableCaching = true;

    // figure out the dbfolder, this is a bit tedious
    const auto calmgr = ExpConfig::Setup::Get(cmd_setup->getValue())->GetCalibrationDataManager();
    const auto calibDataFolder = calmgr->GetCalibrationDataFolder();
    GitInfo gitinfo_db(calibDataFolder);
    DataBase::OnDiskLayout onDiskDB(calibDataFolder);
    if(write && gitinfo_db.IsDirty()) {
        LOG(ERROR) << "Cannot write to dirty database, as undoing will require git checkout";
        return EXIT_FAILURE;
    }

    const auto calibID = cmd_calibration->getValue();
    const auto ranges = [onDiskDB,calibID] () { auto t = onDiskDB.GetDataRanges(calibID); t.sort(); return t; }();
    if(ranges.empty()) {
        cerr << "Could not find any DataRanges for " << calibID << endl;
        return EXIT_FAILURE;
    }

    calibration::gui::AvgBuffer_SavitzkyGolay<TCalibrationData> buffer(cmd_average->getValue(), cmd_sgpol->getValue());

    for(auto& range : ranges) {
        // no need to use the calibration manager,
        // as we already know it's only a DataRange
        // that's different from Ant-calib-dump, for example
        TCalibrationData cdata;
        WrapTFileInput wfi(onDiskDB.GetCurrentFile(range));
        if(wfi.GetObjectClone("cdata",cdata)) {
            if(!check_compatibility(cdata)) {
                LOG(ERROR) << "Incompatible TCalibrationData to previously read found: " << cdata;
                return EXIT_FAILURE;
            }
            // use copy ctor, this is a bit inefficient (but hopefully moving it helps)
            // but it boils down how to handle ownership of ROOT objects (in particular hists)
            // in WrapTFile, then calibration::gui::Manager and then calibration::gui::AvgBuffer
            buffer.Push(make_shared<decltype(cdata)>(std::move(cdata)), range);
        }
        else {
            cerr << "Could not get data for range=" << range << endl;
            return EXIT_FAILURE;
        }
    }

    // indicate no more pushes
    buffer.Flush();

    // before running over the buffer, setup --dump mode
    struct timepoint_t {
        const uint32_t Timestamp;
        const double Value;
        timepoint_t(uint32_t timestamp, double v) :
            Timestamp(timestamp), Value(v) {}
    };

    using timeseries_t = vector<timepoint_t>;
    map<unsigned, timeseries_t> timeseries; // per channel (=key) in map

    if(dump)
        cout << "# Generated with " << std_ext::system::buildCmdLine(argc, argv) << '\n';

    while(!buffer.Empty()) {
        auto& cdata = buffer.CurrentItem();

        if(dump) {
            for(auto& kv : cdata.Data) {
                if(!channels.empty() && !channels.Contains(kv.Key))
                    continue;
                timeseries[kv.Key].emplace_back(cdata.FirstID.Timestamp, kv.Value);
            }
        }
        else if(write) {
            calmgr->Add(cdata, Calibration::AddMode_t::StrictRange);
        }

        buffer.Next();
    }

    if(dump) {
        for(const auto& it_map : timeseries) {
            auto& ch = it_map.first;
            auto& timeseries = it_map.second;
            cout << "# channel=" << ch << '\n';
            for(auto& timepoint : timeseries) {
                cout << timepoint.Timestamp << ' ' << timepoint.Value << '\n';
            }
            cout << '\n' << '\n';
        }
    }
}