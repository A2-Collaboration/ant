#include <map>

#include "detail/taggEffClasses.cc"


#include "base/std_ext/system.h"
#include "base/CmdLine.h"
#include "base/std_ext/string.h"
#include "base/std_ext/time.h"
#include "base/std_ext/math.h"

#include "expconfig/setups/Setup.h"
#include "expconfig/ExpConfig.h"

#include "analysis/plot/root_draw.h"

#include "TTree.h"
#include "TRint.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TMultiGraph.h"

#include "calibration/DataManager.h"
#include "tree/TCalibrationData.h"
#include "calibration/modules/TaggEff.h"

#include "analysis/physics/common/ProcessTaggEff.h"
#include "base/Detector_t.h"



using namespace ant;
using namespace std;
using namespace ant::analysis;
using namespace ant::std_ext;
static volatile bool interrupt = false;

const map<string,TID> startIDs({ {"Setup_2014_07_EPT_Prod", TID(1406592000)},
                                 {"Setup_2014_10_EPT_Prod", TID(1413244800)},
                                 {"Setup_2014_12_EPT_Prod", TID(1417395600)} });

static TH1D* hist_channels(nullptr);
static TH2D* hist_channels_2d(nullptr);
static TH2D* hist_channels_2d_errors(nullptr);

static bool noStore = false;

auto failExit = [] (const string& message)
{
    LOG(ERROR) << message;
    exit(EXIT_FAILURE);
};

void fillHistSingle( const vector<double>& data, const vector<double>& dataErrors);
void fillHistTime( const vector<double>& data, const vector<double>& dataErrors, const unsigned time );

void storeResult(const taggEff_t& result, shared_ptr<calibration::DataManager> manager, const string& calibrationName);

taggEffTriple_t* processFiles(const vector<string>& files);
bool processCSV(const string& csvFile);
void processManualData(const string& dataFile);

int main( int argc, char** argv )
{
    SetupLogger();

    signal(SIGINT, [] (int) {
        LOG(INFO) << ">>> Interrupted";
        interrupt = true;
    });

    TCLAP::CmdLine cmd("Ant-makeTaggEff", ' ', "0.1");

    //modes
    auto mode_csv        = cmd.add<TCLAP::SwitchArg>("","csv",
                                                    "CSV file with bkg1-run-bkg2 groups - calibration data will be added as right open patches for each group");
    auto mode_group      = cmd.add<TCLAP::SwitchArg>("","group",
                                                    "provide single group for adding a right open patch");
    auto mode_manual     = cmd.add<TCLAP::SwitchArg>("","manual",
                                                    "manually set taggeffs per channel for given setup");

    //register modes here:
    auto modes = {&mode_csv, &mode_group, &mode_manual};

    // manual settings:
    auto cmd_setup      = cmd.add<TCLAP::ValueArg<string>>("","setup", "set setup manually",        false, "", "name");
    auto cmd_startDate  = cmd.add<TCLAP::ValueArg<string>>("","start", "set start date manually",   false, "", "date");
    auto cmd_stopDate   = cmd.add<TCLAP::ValueArg<string>>("","stop",  "set end date manually",     false, "", "date");

    //switches
    auto cmd_batchmode  = cmd.add<TCLAP::SwitchArg>("b","batch",  "Run in batch mode (no ROOT shell afterwards)");
    auto cmd_nostore    = cmd.add<TCLAP::SwitchArg>("n","nostore","don't store, only show results");

    //files
    auto cmd_filelist   = cmd.add<TCLAP::UnlabeledMultiArg<string>>("inputfiles","inputfiles to read from",true,"inputfiles");


    cmd.parse(argc, argv);


    auto inputCount = 0u;
    for ( auto m: modes )
        inputCount+=m->get()->isSet();
    if (inputCount!=1)
    {
        string msg = "Exactly one mode is allowed:  ";
        for ( auto m: modes)
            msg += std_ext::formatter() << "--" << m->get()->getName() << "  ";
        failExit(msg);
    }
    auto fileList = cmd_filelist->getValue();
    noStore = cmd_nostore->isSet();

    if (mode_csv->isSet())
    {
        if (fileList.size() != 1)
            failExit("Provide one csv file!");
        processCSV(fileList.at(0));
    }

    taggEffTriple_t* triple_grp = nullptr;
    if (mode_group->isSet())
    {
        if (fileList.size() != 3)
            failExit("Group should contain three files in following order: 1st background, TaggEff-run, background 2.");
        triple_grp = processFiles(fileList);
    }

    if (mode_manual->isSet())
    {
        if (fileList.size() != 1)
            failExit("Provide one data file!");
        processManualData(fileList.at(0));
    }



    argc=1; // prevent TRint to parse any cmdline except prog name
    auto app = cmd_batchmode->isSet() || !std_ext::system::isInteractive()
               ? nullptr
               : std_ext::make_unique<TRint>("Ant-makeSigmas",&argc,argv,nullptr,0,true);

    if(app) {
        canvas c("TaggEff");
        if (hist_channels)
            c << drawoption("E") << hist_channels;
        if (hist_channels_2d)
            c << drawoption("colz") << hist_channels_2d << endr;
        if (hist_channels_2d_errors)
            c << drawoption("colz") << hist_channels_2d_errors << endr;
        c << endc;

        if ( triple_grp )
        {
            canvas control("cotrol");
            auto mg = new TMultiGraph();
            mg->Add(triple_grp->AvgRates);
            mg->Add(triple_grp->avgRatesSub);
            control << drawoption("AP") << mg << endc;
            mg->GetXaxis()->SetTitle("time [s]");
            mg->GetYaxis()->SetTitle("avg. rate [Hz]");

            canvas control_channels("channels");
            control_channels << drawoption("AP");
            for (const auto& b: triple_grp->bkgFits)
                control_channels  << b.Graph;
            control_channels << endc;
            for (const auto& b: triple_grp->bkgFits)
            {
                b.Graph->GetXaxis()->SetTitle("time [s]");
                b.Graph->GetYaxis()->SetTitle("rate [Hz]");
            }
        }

        app->Run(kTRUE); // really important to return...
    }


    return EXIT_SUCCESS;
}

void processManualData(const string& dataFile)
{
    ifstream fstream(dataFile);


    if (!fstream)
        failExit(std_ext::formatter() << "Error opening File " << dataFile << ".");

    auto nChannels = 47u;

    taggEff_t result("setup",TID(),nChannels);

    auto ch = 0u;
    while(fstream)
    {
        string line;
        if (!(getline(fstream,line)))
            break;
        if (ch >= nChannels)
        {
            LOG(WARNING) << "File contains more entries than tagger-channels in setup! skipping the rest";
            break;
        }

        istringstream sstream(line);
        double val(0.);
        sstream >> val;
        result.TaggEffs.emplace_back(val);
        //default error is zero
        val = 0;
        sstream >> val;
        result.TaggEffErrors.emplace_back(val);
        ch++;
    }

    fillHistSingle(result.TaggEffs,result.TaggEffErrors);

    if (!noStore)
    {
        auto manager = ExpConfig::Setup::GetLastFound()->GetCalibrationDataManager();
        storeResult(result,manager,calibration::TaggEff::GetDataName());
    }
}

taggEffTriple_t* processFiles(const vector<string>& files)
{
    auto taggEff = new taggEffTriple_t(files.at(0),files.at(1),files.at(2));
    taggEff_t result = taggEff->GetTaggEffSubtracted();
    auto manager = ExpConfig::Setup::GetLastFound()->GetCalibrationDataManager();

    fillHistSingle(result.TaggEffs,result.TaggEffErrors);

    if (!noStore)
        storeResult(result,manager,calibration::TaggEff::GetDataName());
    return taggEff;
}

bool processCSV(const string& csvFile)
{
    shared_ptr<calibration::DataManager> manager = nullptr;

    ifstream csvStream(csvFile);
    if (!csvStream)
    {
        failExit(std_ext::formatter() << "Error reading File list " << csvFile << ".");
        return false;
    }

    string setupName;
    size_t n_TaggEffs(0);

    while (csvStream && !interrupt )
    {
        string line;
        if (!getline(csvStream, line))
        {
            LOG(INFO) << "Done reading File list " << csvFile << ".";
            break;
        }


        istringstream sstream(line);
        vector<string> record;

        while (sstream)
        {
            string s;
            if (!getline(sstream,s,','))
                break;
            record.emplace_back(s);
        }

        if (record.size() != 3)
            throw runtime_error("Found line with wrong number of files, check your file list.");

        taggEffTriple_t taggEff(record.at(0),record.at(1),record.at(2));
        taggEff_t result = taggEff.GetTaggEffSubtracted();

        //check if setup is valid for this method --> String in map?
        auto it_beamtime = startIDs.find(result.Setup);

        if (it_beamtime == startIDs.end())
            throw runtime_error("Setup not valid for csv mode!");
        if (n_TaggEffs == 0)
        {
            manager = ExpConfig::Setup::GetLastFound()->GetCalibrationDataManager();
            result.FirstID = it_beamtime->second;
        }
        if (n_TaggEffs > 0 && result.Setup != setupName )
                throw runtime_error("Different Setupnames within file list found!");

        fillHistTime(result.TaggEffs,result.TaggEffErrors,result.FirstID.Timestamp);
        if (!noStore)
            storeResult(result,manager,calibration::TaggEff::GetDataName());
        setupName = result.Setup;
        n_TaggEffs++;
    }

    return true;
}



void fillHistSingle( const vector<double>& data, const vector<double>& dataErrors)
{
    if ( !hist_channels )
    {
        hist_channels = new TH1D("taggEff", "Generated Tagging Efficiencies", data.size(),0,data.size());
        hist_channels->SetXTitle("channel");
        hist_channels->SetYTitle("Efficiency");
    }

    for ( auto ch = 0u ; ch < data.size() ; ++ch)
    {
        hist_channels->Fill(ch,data.at(ch));
        hist_channels->SetBinError(ch+1,dataErrors.at(ch)); // bin-numbering ( bin 0 is underflow ) ...
    }
}

void fillHistTime( const vector<double>& data, const vector<double>& dataErrors, const unsigned time )
{
    if ( !hist_channels_2d )
    {
        hist_channels_2d = new TH2D("taggEffTime","Generated Tagging Efficiencies", 1 ,0,0,data.size(),0,data.size());
        hist_channels_2d->SetXTitle("");
        hist_channels_2d->SetYTitle("channel");
    }
    if ( !hist_channels_2d_errors )
    {
        hist_channels_2d_errors = new TH2D("taggEffErrTime","Tagging Efficiencies - relative errors", 1 ,0,0,data.size(),0,data.size());
        hist_channels_2d_errors->SetXTitle("");
        hist_channels_2d_errors->SetYTitle("channel");
    }

    for ( auto ch = 0u ; ch < data.size() ; ++ch )
    {
        hist_channels_2d->Fill(std_ext::to_iso8601(time).c_str(), ch,data.at(ch));
        hist_channels_2d_errors->Fill(std_ext::to_iso8601(time).c_str(), ch,std::abs( 1.0 * dataErrors.at(ch) / data.at(ch)));
    }
}

void storeResult(const taggEff_t& result, shared_ptr<calibration::DataManager> manager, const string& calibrationName)
{
    TCalibrationData cdata(
                calibrationName,
                result.FirstID,
                result.FirstID
                );

    for ( auto ch=0u ; ch < result.TaggEffs.size() ; ++ch )
    {
        cdata.Data.emplace_back(ch,result.TaggEffs.at(ch));
        // see mudule: errors stored in fit-parameters!!!
        cdata.FitParameters.emplace_back(ch,vector<double>(1,result.TaggEffErrors.at(ch)));
    }

    manager->Add(cdata, Calibration::AddMode_t::RightOpen);
}


