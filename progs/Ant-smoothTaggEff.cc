#include "calibration/DataManager.h"
#include "calibration/DataBase.h"
#include "calibration/gui/AvgBuffer.h"

#include "expconfig/ExpConfig.h"

#include "tree/TCalibrationData.h"

#include "tclap/CmdLine.h"
#include "tclap/ValuesConstraintExtra.h"
#include "base/std_ext/system.h"
#include "base/piecewise_interval.h"
#include "base/Logger.h"

#include "detail/tools.h"

using namespace std;
using namespace ant;
using namespace ant::calibration;

using  OnDiskDB_t = DataBase::OnDiskLayout;
using  CalRange_t = OnDiskDB_t::Range_t;

TCalibrationData smoothData(const OnDiskDB_t& onDiskDB, const CalRange_t& range);


int main(int argc, char** argv)
{
    SetupLogger();

    TCLAP::CmdLine cmd("Ant-calib-smooth - applies low-pass filter on calibration data ranges", ' ', "0.1");

    TCLAP::ValuesConstraintExtra<decltype(ExpConfig::Setup::GetNames())> allowedsetupnames(ExpConfig::Setup::GetNames());
    auto cmd_setup  = cmd.add<TCLAP::ValueArg<string>>("s","setup","Use setup to determine calibration database path",true,"", &allowedsetupnames);

    auto cmd_polOrd = cmd.add<TCLAP::ValueArg<unsigned>>("","polyorder","Polynom order for Savitzky-Golay filter (zero is moving average)", false, 4, "polorder");

    auto cmd_dump  = cmd.add<TCLAP::SwitchArg>("","dump","Dump to stdout",false);

    auto cmd_write  = cmd.add<TCLAP::SwitchArg>("","write","Write to database",false);

    cmd.parse(argc, argv);

    if(!(cmd_dump->isSet() ^ cmd_write->isSet())) {
        LOG(ERROR) << "Use either --dump or --write, but not both.";
        return EXIT_FAILURE;
    }

//    const auto dump = cmd_dump->getValue();
    const auto write = cmd_write->getValue();

    // enable caching of the calibration database
    DataBase::OnDiskLayout::EnableCaching = true;

    // figure out the dbfolder, this is a bit tedious
    ExpConfig::Setup::SetByName(cmd_setup->getValue());
    const auto calmgr = ExpConfig::Setup::Get().GetCalibrationDataManager();
    const auto calibDataFolder = calmgr->GetCalibrationDataFolder();
    GitInfo gitinfo_db(calibDataFolder);
    DataBase::OnDiskLayout onDiskDB(calibDataFolder);
    if(write && gitinfo_db.IsDirty()) {
        LOG(ERROR) << "Cannot write to dirty database, as undoing will require git checkout";
        return EXIT_FAILURE;
    }

    const auto calibID = "TaggEff";
    const auto ranges = [onDiskDB,calibID] () { auto t = onDiskDB.GetDataRanges(calibID); t.sort(); return t; }();
    if(ranges.empty()) {
        cerr << "Could not find any DataRanges for " << calibID << endl;
        return EXIT_FAILURE;
    }


    for(auto& range : ranges) {
        const auto smoothedData = smoothData(onDiskDB,range);

    }


}

TCalibrationData smoothData(const OnDiskDB_t& onDiskDB, const CalRange_t& range)
{
    TCalibrationData cdata;
    WrapTFileInput wfi(onDiskDB.GetCurrentFile(range));
    if(!wfi.GetObjectClone("cdata",cdata))
    {
        string errmsg = std_ext::formatter() << "Could not get data for range=" << range;
        throw runtime_error(errmsg);
    }
    return cdata;
}


