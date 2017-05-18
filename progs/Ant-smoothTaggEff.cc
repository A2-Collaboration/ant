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
#include "base/PlotExt.h"

#include "detail/tools.h"

#include "TGraphErrors.h"
#include "TF1.h"

#include <iostream>

using namespace std;
using namespace ant;
using namespace ant::calibration;

using  OnDiskDB_t = DataBase::OnDiskLayout;
using  CalRange_t = OnDiskDB_t::Range_t;

TCalibrationData GetData(const OnDiskDB_t& onDiskDB, const CalRange_t& range);
TGraphErrors* toCalGraph(const TCalibrationData& cdata);

struct taggEffData_t{
        size_t Channel;
        double TaggEff;
        double TaggEffError;
        taggEffData_t(size_t ch, double taggEff, double taggEffError):
            Channel(ch), TaggEff(taggEff), TaggEffError(taggEffError){}
};

template<typename F>
void forEachChannel(const TCalibrationData& cdata, F func)
{
    for (const auto& data: cdata.Data)
    {
        taggEffData_t ted(data.Key, data.Value, std_ext::NaN);

        for ( const auto& fitP: cdata.FitParameters)
        {
            if (ted.Channel == fitP.Key)
            {
                ted.TaggEffError = fitP.Value.front();
                break;
            }
        }
        func(ted);
    }
}

int main(int argc, char** argv)
{
    SetupLogger();

    TCLAP::CmdLine cmd("Ant-calib-smooth - applies low-pass filter on calibration data ranges", ' ', "0.1");

    TCLAP::ValuesConstraintExtra<decltype(ExpConfig::Setup::GetNames())> allowedsetupnames(ExpConfig::Setup::GetNames());
    auto cmd_setup  = cmd.add<TCLAP::ValueArg<string>>("s","setup","Use setup to determine calibration database path",true,"", &allowedsetupnames);

    auto cmd_calibName  = cmd.add<TCLAP::ValueArg<string>>("c","calibration","Name of calibration, taggeff is stored in",true,"","calibration");

    auto cmd_polOrd = cmd.add<TCLAP::ValueArg<unsigned>>("","polyorder","Polynom order for Savitzky-Golay filter (zero is moving average)", false, 4, "polorder");

    auto cmd_dump  = cmd.add<TCLAP::SwitchArg>("","dump","Dump to stdout",false);

    auto cmd_write  = cmd.add<TCLAP::SwitchArg>("","write","Write to database",false);

    cmd.parse(argc, argv);

    const auto   dump   = cmd_dump->getValue();
    const auto   write  = cmd_write->getValue();
    const string fitstr = std_ext::formatter() << "pol" <<cmd_polOrd->getValue();

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

    const auto calibID = cmd_calibName->getValue();
    const auto ranges = [onDiskDB,calibID] () { auto t = onDiskDB.GetDataRanges(calibID); t.sort(); return t; }();
    if(ranges.empty()) {
        cerr << "Could not find any DataRanges for " << calibID << endl;
        return EXIT_FAILURE;
    }

    for(auto& range : ranges)
    {
        auto cData = GetData(onDiskDB,range);
        auto graph = toCalGraph(cData);
        TF1* fkt = new TF1("fit", fitstr.c_str(), 0, graph->GetN());
        graph->Fit(fkt,"Q");
        string msg = std_ext::formatter() << "Range: " << range << "\n";
        for (auto& data: cData.Data)
        {
            msg += std_ext::formatter() << data.Key   << " "
                                        << data.Value << " ";
            data.Value = fkt->Eval(data.Key);
            msg += std_ext::formatter() << data.Value << "\n";
        }
        if (dump)
            cout << msg << endl;
        if (write)
            calmgr->Add(cData, Calibration::AddMode_t::RightOpen);  //TaggEffs are all right open if not default

    }




}



TCalibrationData GetData(const OnDiskDB_t& onDiskDB, const CalRange_t& range)
{
    TCalibrationData cdata;
    WrapTFileInput wfi(onDiskDB.GetCurrentFile(range));
    if(!wfi.GetObjectClone("cdata",cdata))
            throw runtime_error(std_ext::formatter() << "Could not get data for range = " << range);
    return cdata;
}

TGraphErrors* toCalGraph(const TCalibrationData& cdata)
{
    auto g = new TGraphErrors();
    forEachChannel(
                cdata,
                [g](const taggEffData_t& ted)
                {
                    GraphExt::FillGraphErrors(g,ted.Channel,ted.TaggEff,0,ted.TaggEffError);
                });
    return g;
}



