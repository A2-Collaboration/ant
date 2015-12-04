#include <iostream>
#include <string>
#include <sstream>
#include <vector>


#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TRint.h"
#include "TStyle.h"

#include "analysis/plot/root_draw.h"

#include "base/CmdLine.h"
#include "base/std_ext/string.h"
#include "base/Logger.h"
#include "base/detail/tclap/ValuesConstraintExtra.h"
#include "base/WrapTFile.h"

#include "calibration/DataBase.h"
#include "calibration/DataManager.h"



using namespace std;
using namespace ant;

void print_help();
void list_calibrations(string& cfolder);
TH2D* show_calibration(const string& cfolder, const string& calibID);


int main( int argc, char** argv )
{
    SetupLogger();

    if (argc <= 2)
    {
        print_help();
        return 1;
    }

    string todo(argv[1]);
    string cfolder(argv[2]);

    if ( todo.compare("list") == 0 )
    {
        list_calibrations(cfolder);
        return 0;
    }
    if ( todo.compare("show") == 0 )
    {
//        if (argc != 5)
        if (argc != 4)
        {
            print_help();
            return 1;
        }
        TH2D* plot = show_calibration(cfolder,string(argv[3]));
        if (plot)
        {
            int fake_argc=1;
            char* fake_argv[2];
            fake_argv[0] = argv[0];
            auto app = new TRint("Ant-calib-viewer",&fake_argc,fake_argv);
            gStyle->SetOptStat(kFALSE);

            canvas c("view");
            c << drawoption("colz") << plot << endc;
            plot->SetShowProjectionX(1);


            app->Run(kTRUE);

            return 0;
        }
        cout << "  No Ranges defined for this calibration!" << endl;
    }

    print_help();
    return 1;

}


void print_help()
{
    cout << " Usage:  Ant-calib-viewer list < calibration folder >" << endl
         << "                              --  List calibrationIDs in database" << endl << endl
    //todo visualize MC and defaults...
//         << "                          show < calibration folder > < calibrationID > { default | mc | ranges}" << endl
         << "                          show < calibration folder > < calibrationID >" << endl
         << "                              --  Plot Calibration" << endl << endl ;
}

void list_calibrations(std::string& cfolder)
{
    calibration::DataManager dm(cfolder);
    for (const auto& cID: dm.GetCalibrationIDs())
        cout << cID << endl;
}

void GetData(const calibration::DataBase::OnDiskLayout& onDiskDB, const calibration::DataBase::OnDiskLayout::Range_t& range, TCalibrationData& dataBuffer )
{
    WrapTFileInput wfi(onDiskDB.GetCurrentFile(range));
    wfi.GetObjectClone("cdata",dataBuffer);
}

TH2D* show_calibration(const string& cfolder, const string& calibID)
{

    TH2D* hist = nullptr;

    TCalibrationData dataBuffer;
    calibration::DataBase::OnDiskLayout onDiskDB(cfolder);


    //todo visualize MC and defaults...
    /*
    if (type.compare("default") == 0 )
    {

    }
    if (type.compare("mc") == 0 )
    {
    }
    if (type.compare("ranges") == 0 )
    */

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

    return hist;
}
