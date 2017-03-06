#include "calibration/DataManager.h"
#include "calibration/DataBase.h"

#include "expconfig/ExpConfig.h"

#include "tree/TCalibrationData.h"

#include "base/CmdLine.h"
#include "base/detail/tclap/ValuesConstraintExtra.h"
#include "base/std_ext/system.h"
#include "base/piecewise_interval.h"
#include "base/Logger.h"


using namespace std;
using namespace ant;
using namespace ant::calibration;

TCalibrationData read_cdata(std::istream& is) {
    TCalibrationData cdata;
    string line;
    unsigned n_line = 0;
    while(getline(is, line)) {
        n_line++;
        if(std_ext::string_starts_with(line, "#"))
            continue;
        stringstream ss(line);
        unsigned ch;
        if(!(ss >> ch)) {
            LOG(WARNING) << "Cannot parse channel in line "
                         << n_line << ": " << line;
            continue;
        }
        auto it_item = find_if(cdata.Data.begin(), cdata.Data.end(), [ch] (const TKeyValue<double>& kv) {
            return kv.Key == ch;
        });
        if(it_item != cdata.Data.end()) {
            LOG(WARNING) << "Channel " << ch << " already seen in line "
                         << n_line << ": " << line;
            continue;
        }

        std::vector<double> values;
        {
            double val;
            while(ss >> val)
                values.push_back(val);
        }
        if(values.empty()) {
            LOG(WARNING) << "No values found in line "
                         << n_line << ": " << line;
            continue;
        }

        // add first value as data, the rest as fit parameters, if any (similar to Ant-calib-dump)
        cdata.Data.emplace_back(ch, values.front());
        values.erase(values.begin());
        if(!values.empty())
            cdata.FitParameters.emplace_back(ch, move(values));
    }
    return cdata;
}

int main(int argc, char** argv)
{
    SetupLogger();

    TCLAP::CmdLine cmd("Ant-calib-readin - add text data to calibration database", ' ', "0.1");

    auto cmd_inputfile  = cmd.add<TCLAP::ValueArg<string>>("i","input","Input filename (STDIN by default)",false,"-","inputfile");
    TCLAP::ValuesConstraintExtra<decltype(ExpConfig::Setup::GetNames())> allowedsetupnames(ExpConfig::Setup::GetNames());
    auto cmd_setup  = cmd.add<TCLAP::ValueArg<string>>("s","setup","Use setup to determine calibration database path",true,"", &allowedsetupnames);
    auto cmd_calibration = cmd.add<TCLAP::ValueArg<string>>("c","calibration","Calibration ID", true, "","calibration");
    auto cmd_type_mc  = cmd.add<TCLAP::SwitchArg>("","mc","Add as MC values",false);

    cmd.parse(argc, argv);

    // read the file
    TCalibrationData cdata;
    const auto& inputfile = cmd_inputfile->getValue();
    if(inputfile != "-") {
        ifstream is(inputfile);
        if(!is) {
            LOG(ERROR) << "Cannot open inputfile " << inputfile;
            return EXIT_FAILURE;
        }
        cdata = read_cdata(is);
    }
    else {
        cdata = read_cdata(std::cin);
    }

    if(cdata.Data.empty()) {
        LOG(ERROR) << "Nothing read in";
        return EXIT_FAILURE;
    }

    // set some necessary cdata fields
    cdata.CalibrationID = cmd_calibration->getValue();

    GitInfo info;
    cdata.Author = info.GetUser();

    TID fake_tid(0,0,
                 cmd_type_mc->isSet() ?
                     list<TID::Flags_t>{TID::Flags_t::AdHoc, TID::Flags_t::MC } :
                     list<TID::Flags_t>{TID::Flags_t::AdHoc});
    cdata.FirstID = fake_tid;
    cdata.LastID  = fake_tid;

    // figure out the calib manager and try to add the data as default
    const auto calmgr = ExpConfig::Setup::Get(cmd_setup->getValue())->GetCalibrationDataManager();
    calmgr->Add(cdata, Calibration::AddMode_t::AsDefault);

    return EXIT_SUCCESS;
}