#include "calibration/DataManager.h"
#include "calibration/DataBase.h"

#include "expconfig/ExpConfig.h"

#include "tree/TCalibrationData.h"

#include "tclap/CmdLine.h"
#include "tclap/ValuesConstraintExtra.h"
#include "base/std_ext/system.h"
#include "base/std_ext/container.h"
#include "base/piecewise_interval.h"
#include "base/Logger.h"

#include "detail/tools.h"

using namespace std;
using namespace ant;
using namespace ant::calibration;

int main(int argc, char** argv)
{
    SetupLogger();
    // make logger silent, as Database uses LOG(INFO) to output successful load message
    // which destroys the parseable output on std::cout
    el::Configurations loggerConf;
    loggerConf.setGlobally(el::ConfigurationType::Enabled, "false");
    el::Loggers::reconfigureLogger("default", loggerConf);


    TCLAP::CmdLine cmd("Ant-calib-dump - dump data from calibration database as text", ' ', "0.1");

    TCLAP::ValuesConstraintExtra<decltype(ExpConfig::Setup::GetNames())> allowedsetupnames(ExpConfig::Setup::GetNames());
    auto cmd_setup  = cmd.add<TCLAP::ValueArg<string>>("s","setup","Use setup to determine calibration database path",true,"", &allowedsetupnames);
    auto cmd_calibration = cmd.add<TCLAP::ValueArg<string>>("c","calibration","Calibration ID", true, "","calibration");
    auto cmd_type_mc  = cmd.add<TCLAP::SwitchArg>("","mc","Dump MC values (time-independent)",false);
    auto cmd_type_datadefault  = cmd.add<TCLAP::SwitchArg>("","datadefault","Dump DataDefault values (time-independent)",false);
    auto cmd_channels = cmd.add<TCLAP::MultiArg<string>>("","ch","Ranges of channels, e.g. 400-412",false,"channels");
    auto cmd_params = cmd.add<TCLAP::MultiArg<string>>("","params","Index for value=0/FitParams=1,2,... as columns ",false,"params");
    auto cmd_noignore  = cmd.add<TCLAP::SwitchArg>("","noignore","Do not output ignored channels",false);
    auto cmd_notouches  = cmd.add<TCLAP::SwitchArg>("","notouches","Do not output channels which touch hole",false);

    cmd.parse(argc, argv);

    // do some options handling
    if(cmd_type_mc->isSet() && cmd_type_datadefault->isSet()) {
        cerr << "Flags --mc and --datadefault are mutually exclusive" << endl;
        return EXIT_FAILURE;
    }
    using datatype_t = DataBase::OnDiskLayout::Type_t;
    datatype_t dataType = datatype_t::DataRanges;
    if(cmd_type_mc->isSet())
        dataType = datatype_t::MC;
    if(cmd_type_datadefault->isSet())
        dataType = datatype_t::DataDefault;

    const auto channels = progs::tools::parse_cmdline_ranges(cmd_channels->getValue());
    const auto params = progs::tools::parse_cmdline_ranges(cmd_params->getValue());

    const bool noIgnored = cmd_noignore->isSet();
    const bool noTouchesHole = cmd_notouches->isSet();

    // enable caching of the calibration database
    DataBase::OnDiskLayout::EnableCaching = true;

    // figure out the dbfolder

    ExpConfig::Setup::SetByName(cmd_setup->getValue());
    const auto calmgr = ExpConfig::Setup::Get().GetCalibrationDataManager();
    DataBase::OnDiskLayout onDiskDB(calmgr->GetCalibrationDataFolder());

    const auto calibID = cmd_calibration->getValue();
    auto available_calibrations = calmgr->GetCalibrationIDs();
    if(!std_ext::contains(available_calibrations, calibID)) {
        cerr << "Cannot find calibration ID " << calibID << endl;
        return EXIT_FAILURE;
    }
    const auto dummyTID = dataType == datatype_t::MC ? TID(0,0,{TID::Flags_t::MC}) : TID(0,0);
    const auto ranges = dataType == datatype_t::DataRanges ?
                            [onDiskDB,calibID] () { auto t = onDiskDB.GetDataRanges(calibID); t.sort(); return t; }() :
                            list<DataBase::OnDiskLayout::Range_t>{{{dummyTID, dummyTID},""}}; // one dummy element, but MC flag is important
    if(ranges.empty()) {
        cerr << "Could not find ranges" << endl;
        return EXIT_FAILURE;
    }

    // ugly way of getting the detector pointer, we infer it from the CalibId,
    // which might fail in general
    auto detector_string = std_ext::tokenize_string(calibID,"_").front();
    auto det = ExpConfig::Setup::GetDetector(Detector_t::FromString(detector_string));
    auto cldet = dynamic_pointer_cast<ClusterDetector_t>(det);
    if(noTouchesHole && cldet==nullptr) {
        cerr << "You can only use --notouches for cluster detectors" << endl;
        return EXIT_FAILURE;
    }

    const auto is_channel_excluded = [channels, det, cldet, noIgnored, noTouchesHole] (const unsigned ch) {
        if(!channels.empty() && !channels.Contains(ch))
            return true;
        if(noIgnored && det->IsIgnored(ch))
            return true;
        if(noTouchesHole && cldet->GetClusterElement(ch)->TouchesHole)
            return true;
        return false;
    };

    struct timepoint_t {
        uint32_t Timestamp;
        vector<double> ValueFitParams;
        timepoint_t(uint32_t timestamp, const vector<double>& p) :
            Timestamp(timestamp), ValueFitParams(p) {}
    };

    using timeseries_t = vector<timepoint_t>;
    map<unsigned, timeseries_t> timeseries; // per channel (=key) in map

    for(auto& range : ranges) {
        TCalibrationData cdata;
        if(calmgr->GetData(calibID, range.Start(), cdata)) {
            // data first
            if(params.empty() || params.Contains(0)) {
                for(auto& kv : cdata.Data) {
                    if(is_channel_excluded(kv.Key))
                        continue;
                    timeseries[kv.Key].emplace_back(range.Start().Timestamp, vector<double>{kv.Value});
                }
            }

            // then fit params
            for(auto& kv : cdata.FitParameters) {
                if(is_channel_excluded(kv.Key))
                    continue;
                // need to check if data has already put next item in timeseries,
                // then just extend that
                timeseries_t& series = timeseries[kv.Key];
                if(series.empty() || series.back().Timestamp != range.Start().Timestamp)
                    series.emplace_back(range.Start().Timestamp, vector<double>{});
                auto& valueFitParams = series.back().ValueFitParams;

                for(int p=0;p<int(kv.Value.size());p++) {
                    if(params.empty() || params.Contains(p+1)) // offset by one as Data is index=0
                        valueFitParams.push_back(kv.Value[p]);
                }
            }

        }
        else {
            cerr << "Could not get data for range=" << range << endl;
            return EXIT_FAILURE;
        }
    }

    cout << "# Generated with " << std_ext::system::buildCmdLine(argc, argv) << '\n';

    if(dataType == datatype_t::DataRanges) {
        // each channel has a timeseries of data
        // so time is x coordinate for plotting
        for(const auto& it_map : timeseries) {
            auto& ch = it_map.first;
            auto& timeseries = it_map.second;
            cout << "# channel=" << ch << '\n';

            for(auto& timepoint : timeseries) {
                cout << timepoint.Timestamp;
                for(auto& v : timepoint.ValueFitParams)
                    cout << ' ' << v;
                cout << '\n';
            }

            cout << '\n' << '\n';
        }
    }
    else {
        // default/MC does not have a timeseries,
        // so dump channel as x coordinate
        for(const auto& it_map : timeseries) {
            auto& ch = it_map.first;
            auto& timeseries = it_map.second;
            cout << ch;
            for(auto& v : timeseries.front().ValueFitParams)
                cout << ' ' << v;
            cout << '\n';
        }
    }
}