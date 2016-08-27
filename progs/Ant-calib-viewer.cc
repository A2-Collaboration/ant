#include <iostream>
#include <string>
#include <sstream>
#include <vector>



#include "analysis/plot/root_draw.h"

#include "base/CmdLine.h"
#include "base/std_ext/string.h"
#include "base/std_ext/system.h"
#include "base/Logger.h"
#include "base/detail/tclap/ValuesConstraintExtra.h"
#include "base/WrapTFile.h"

#include "calibration/DataBase.h"
#include "calibration/DataManager.h"

#include "expconfig/ExpConfig.h"

#include "tree/TCalibrationData.h"

#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TRint.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TTree.h"


using namespace std;
using namespace ant;
using namespace ant::calibration;

void show_2d(const string& dbfolder, const string& calibID);
void show_time(const string& dbfolder, const string& calibID);


int main(int argc, char** argv)
{
    SetupLogger();

    TCLAP::CmdLine cmd("Ant-calib-viewer - plot calibration parameters from database", ' ', "0.1");
    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"level");

    TCLAP::ValuesConstraintExtra<decltype(ExpConfig::Setup::GetNames())> allowedsetupnames(ExpConfig::Setup::GetNames());
    auto cmd_setup  = cmd.add<TCLAP::ValueArg<string>>("s","setup","Use setup to determine calibration database path",false,"", &allowedsetupnames);
    auto cmd_dbfolder = cmd.add<TCLAP::ValueArg<string>>("d","dbfolder","Path to calibration database, with trailing 'calibration'", false, "", "dbfolder");
    auto cmd_calibration = cmd.add<TCLAP::ValueArg<string>>("c","calibration","Calibration ID", true, "","calibration");

    TCLAP::SwitchArg cmd_mode_show_2d("","show_2d","Show 2D Hist with SetProjectionX", false);
    TCLAP::SwitchArg cmd_mode_show_convergence("","show_convergence","Show values vs. iteration", false);
    TCLAP::SwitchArg cmd_mode_show_time("","show_time","Show values vs. time", false);

    cmd.xorAdd({&cmd_mode_show_2d, &cmd_mode_show_convergence, &cmd_mode_show_time});

    cmd.parse(argc, argv);

    if(cmd_verbose->isSet())
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());

    // figure out the dbfolder
    string dbfolder;
    if(cmd_dbfolder->isSet()) {
        dbfolder = cmd_dbfolder->getValue();
    }
    else if (cmd_setup->isSet()) {
        ExpConfig::Setup::SetManualName(cmd_setup->getValue());
        auto calmgr = ExpConfig::Setup::GetLastFound()-> GetCalibrationDataManager();
        dbfolder = calmgr->GetCalibrationDataFolder();
    }
    else {
        LOG(ERROR) << "Neither setupname nor dbfolder specified.";
        return 1;
    }

    // check if calibration ID exists at least
    string calibrationID = cmd_calibration->getValue();
    auto folders = std_ext::system::lsFiles(dbfolder+"/"+calibrationID,"",true,false);
    if(folders.empty()) {
        LOG(ERROR) << "Calibration ID '" << calibrationID << "' does not exist in " << dbfolder;
        return 1;
    }

    int fake_argc=0;
    char* fake_argv[1];
    fake_argv[0] = argv[0];
    auto app = new TRint("Ant-calib-viewer",&fake_argc,fake_argv,nullptr,0,true);

    if(cmd_mode_show_2d.isSet()) {
        show_2d(dbfolder, calibrationID);
    }

    if(cmd_mode_show_convergence.isSet()) {
        LOG(INFO) << "Not implemented yet";
    }

    if(cmd_mode_show_time.isSet()) {
        show_time(dbfolder, calibrationID);
    }

    app->Run(kTRUE);

    return 0;
}

void GetData(const DataBase::OnDiskLayout& onDiskDB,
             const DataBase::OnDiskLayout::Range_t& range,
             TCalibrationData& dataBuffer )
{
    WrapTFileInput wfi(onDiskDB.GetCurrentFile(range));
    wfi.GetObjectClone("cdata",dataBuffer);
}

void show_time(const string& dbfolder, const string& calibID)
{
    TCalibrationData dataBuffer;
    DataBase::OnDiskLayout onDiskDB(dbfolder);


    auto points_x = new std::vector<double>();
    auto points_y = new std::vector<double>();

    auto drawTree = new TTree();
    double key;
    double timestamp;
    double value;
    drawTree->Branch("key",addressof(key));
    drawTree->Branch("time",addressof(timestamp));
    drawTree->Branch("value",addressof(value));


    auto ranges = onDiskDB.GetDataRanges(calibID);

    unsigned i = 0;
    for (const auto& range: ranges)
    {
        GetData(onDiskDB, range, dataBuffer);


        const double length = dataBuffer.Data.size();
        const double span = 0.7;
        unsigned j = 0;
        for(const TCalibrationData::Entry& entry: dataBuffer.Data) {
            const double x = i + span*(double)j/length - span/2;
            points_x->push_back(x);
            points_y->push_back(entry.Value);

            key   = (double) j;
            timestamp = 1.0 * ((double)dataBuffer.FirstID.Timestamp + (double)dataBuffer.LastID.Timestamp) / 2.0;
            value  = entry.Value;
            drawTree->Fill();

            j++;
        }

        i++;
    }

    auto graph = new TGraph(points_x->size(), points_x->data(), points_y->data());
    graph->SetTitle(calibID.c_str());
    graph->GetXaxis()->SetTitle("Range");
    graph->GetYaxis()->SetTitle("Value");


    canvas c("view");
    c << drawoption("AP") << graph << drawoption("PCOL") << TTree_drawable(drawTree,"value:time:k") << endc;

}

void show_2d(const string& dbfolder, const string& calibID)
{

    TH2D* hist = nullptr;

    TCalibrationData dataBuffer;
    DataBase::OnDiskLayout onDiskDB(dbfolder);


    //todo visualize MC and defaults...

    auto ranges = onDiskDB.GetDataRanges(calibID);

    unsigned i = 0;
    for (const auto& range: ranges)
    {
        GetData(onDiskDB, range, dataBuffer);

        if (!hist)
        {
            hist = new TH2D("ranges","",ranges.size(),0,ranges.size(), dataBuffer.Data.size(),0,dataBuffer.Data.size());
            hist->SetXTitle("ranges");
            hist->SetYTitle("channel");
        }

        for (const auto& entry: dataBuffer.Data)
            hist->Fill(i,entry.Key,entry.Value);
        i++;
    }
    if(!hist)
        return;


    gStyle->SetOptStat(kFALSE);

    canvas c("view");
    c << drawoption("colz") << hist << endc;
    hist->SetShowProjectionX(1);


}
