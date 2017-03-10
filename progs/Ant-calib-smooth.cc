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

int main(int argc, char** argv)
{
    SetupLogger();

    TCLAP::CmdLine cmd("Ant-calib-smooth - applies low-pass filter on calibration data ranges", ' ', "0.1");

    TCLAP::ValuesConstraintExtra<decltype(ExpConfig::Setup::GetNames())> allowedsetupnames(ExpConfig::Setup::GetNames());
    auto cmd_setup  = cmd.add<TCLAP::ValueArg<string>>("s","setup","Use setup to determine calibration database path",true,"", &allowedsetupnames);
    auto cmd_calibration = cmd.add<TCLAP::ValueArg<string>>("c","calibration","Calibration ID", true, "","calibration");

    cmd.parse(argc, argv);

    // enable caching of the calibration database
    DataBase::OnDiskLayout::EnableCaching = true;

    // figure out the dbfolder

    const auto calmgr = ExpConfig::Setup::Get(cmd_setup->getValue())->GetCalibrationDataManager();
    DataBase::OnDiskLayout onDiskDB(calmgr->GetCalibrationDataFolder());

    const auto calibID = cmd_calibration->getValue();
    const auto ranges = [onDiskDB,calibID] () { auto t = onDiskDB.GetDataRanges(calibID); t.sort(); return t; }();
    if(ranges.empty()) {
        cerr << "Could not find DataRanges for " << calibID << endl;
        return EXIT_FAILURE;
    }

    for(auto& range : ranges) {
        TCalibrationData cdata;
        if(calmgr->GetData(calibID, range.Start(), cdata)) {


        }
        else {
            cerr << "Could not get data for range=" << range << endl;
            return EXIT_FAILURE;
        }
    }
}