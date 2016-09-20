#include <string>
#include <map>

#include "analysis/physics/common/ProcessTaggEff.h"
#include "analysis/plot/root_draw.h"

#include "base/CmdLine.h"
#include "base/Detector_t.h"
#include "base/Logger.h"
#include "base/std_ext/math.h"
#include "base/std_ext/string.h"
#include "base/std_ext/system.h"
#include "base/std_ext/time.h"

#include "calibration/DataManager.h"
#include "calibration/modules/TaggEff.h"

#include "detail/taggEffClasses.h"

#include "expconfig/ExpConfig.h"
#include "expconfig/setups/Setup.h"

#include "tree/TCalibrationData.h"


#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMultiGraph.h"
#include "TRint.h"
#include "TTree.h"
#include "TGraphErrors.h"



using namespace ant;
using namespace std;
using namespace ant::analysis;
using namespace ant::std_ext;
using namespace ant::progs::taggeff;




const map<string,TID> startIDs({ {"Setup_2014_07_EPT_Prod", TID(1406592000)},
                                 {"Setup_2014_10_EPT_Prod", TID(1413244800)},
                                 {"Setup_2014_12_EPT_Prod", TID(1417395600)} });

constexpr double chi2cut_channels = 15.0;

static TH1D* hist_channels(nullptr);
static TH2D* hist_channels_2d(nullptr);
static TH2D* hist_channels_2d_errors(nullptr);
static TGraphErrors* graphTime(nullptr);

static shared_ptr<HistogramFactory> histfac = nullptr;

static volatile bool interrupt = false;
static bool noStore = false;
static bool histOut = false;

auto failExit = [] (const string& message)
{
    LOG(ERROR) << message;
    exit(EXIT_FAILURE);
};

void fillHistSingle( const vector<double>& data, const vector<double>& dataErrors);
void fillHistTime( const vector<double>& data, const vector<double>& dataErrors, const unsigned time );

void storeHist(const taggEff_t& result);
void storeResult(const taggEff_t& result, shared_ptr<calibration::DataManager> manager, const string& calibrationName,
                 const Calibration::AddMode_t& addMode = Calibration::AddMode_t::RightOpen);

taggEffTriple_t* processFiles(const vector<string>& files);
void mediateCSV(const vector<string>& csvFiles);
void processCSV(const string& csvFile);
void processManualData(const string& dataFile, const string& setupName);

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
    auto mode_csv_mean   = cmd.add<TCLAP::SwitchArg>("","csv-mean",
                                                     "CSV file with bkg1-run-bkg2 groups - calibration data will be the mean over all measurements.");
    auto mode_group      = cmd.add<TCLAP::SwitchArg>("","group",
                                                    "provide single group for adding a right open patch");
    auto mode_manual     = cmd.add<TCLAP::SwitchArg>("","manual",
                                                    "manually set taggeffs per channel for given setup");

    //register modes here:
    auto modes = {&mode_csv, &mode_csv_mean, &mode_group, &mode_manual};

    // other settings:
    auto cmd_setup      = cmd.add<TCLAP::ValueArg<string>>("","setup", "set setup manually",        false, "", "name");

    //switches
    auto cmd_batchmode  = cmd.add<TCLAP::SwitchArg>("b", "batch",     "Run in batch mode (no ROOT shell afterwards)");
    auto cmd_nostore    = cmd.add<TCLAP::SwitchArg>("n", "nostore",   "don't store, only show results");
    auto cmd_histout    = cmd.add<TCLAP::SwitchArg>("",  "save-hist", "Save results to a histogram");

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
    histOut = cmd_histout->isSet();

    if (mode_csv->isSet())
    {
        if (cmd_histout)
        if (fileList.size() != 1)
            failExit("Provide one csv file!");
        processCSV(fileList.at(0));
    }

    if (mode_csv_mean->isSet())
    {
        if (cmd_histout)
        if (fileList.size() < 1)
            failExit("Provide at least one csv file!");
        mediateCSV(fileList);
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
        if (!cmd_setup->isSet())
            failExit("Must provide setup name for data.");
        if (fileList.size() != 1)
            failExit("Provide one data file!");
        processManualData(fileList.at(0),cmd_setup->getValue());
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

        if ( graphTime )
        {
            canvas("means") << graphTime << endc;
            graphTime->GetXaxis()->SetTimeDisplay(1);
            graphTime->GetXaxis()->SetTimeFormat("%d.%m.%Y %F 1970-01-01 00:00:00");
        }

        if ( triple_grp )
        {
            canvas control("cotrol");
            auto mg = new TMultiGraph();
            mg->Add(triple_grp->AvgBkgRates);
            mg->Add(triple_grp->avgRatesSub);
            mg->Add(triple_grp->avgRates);
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

void processManualData(const string& dataFile, const string& setupName)
{
    ifstream fstream(dataFile);

    if (!fstream)
        failExit(std_ext::formatter() << "Error opening File " << dataFile << ".");

    ExpConfig::Setup::SetManualName(setupName);
    auto nChannels = ExpConfig::Setup::GetLastFound()->GetDetector<TaggerDetector_t>()->GetNChannels();
    taggEff_t result("setup",TID(),nChannels);

    auto ch = 0u;
    while(fstream)
    {
        string line;
        if (!(getline(fstream,line)))
            break;
        if (ch >= nChannels)
            failExit("File contains more entries than tagger-channels in setup!");

        istringstream sstream(line);
        double val(0.);
        sstream >> val;
        result.TaggEffs.at(ch) = val;
        //default error is zero
        val = 0;
        sstream >> val;
        result.TaggEffErrors.at(ch) = val;
        ch++;
    }

    if (ch != nChannels)
        failExit("File contains less entries than tagger-channels in setup.");

    fillHistSingle(result.TaggEffs,result.TaggEffErrors);

    if (!noStore)
    {
        auto manager = ExpConfig::Setup::GetLastFound()->GetCalibrationDataManager();
        storeResult(result,manager,calibration::TaggEff::GetDataName());
    }
    if ( histOut )
        storeHist(result);
}

taggEffTriple_t* processFiles(const vector<string>& files)
{
    auto taggEff = new taggEffTriple_t(files.at(0),files.at(1),files.at(2));
    taggEff_t result = taggEff->GetTaggEffSubtracted();
    auto manager = ExpConfig::Setup::GetLastFound()->GetCalibrationDataManager();

    fillHistSingle(result.TaggEffs,result.TaggEffErrors);

    if (!noStore)
        storeResult(result,manager,calibration::TaggEff::GetDataName());
    if (histOut)
        storeHist(result);

    return taggEff;
}

void processCSV(const string& csvFile)
{
    shared_ptr<calibration::DataManager> manager = nullptr;

    ifstream csvStream(csvFile);
    if (!csvStream)
        failExit(std_ext::formatter() << "Error reading File list " << csvFile << ".");


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
}

void mediateCSV(const vector<string>& csvFiles)
{
    string setupName;
    size_t n_TaggEffs(0);
    size_t nCh(0);

    shared_ptr<calibration::DataManager> manager = nullptr;

    vector<double> vS;
    vector<double> vSy;

    for ( const auto& csvFile: csvFiles)
    {
        ifstream csvStream(csvFile);
        if (!csvStream)
            failExit(std_ext::formatter() << "Error reading File list " << csvFile << ".");

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
            if (n_TaggEffs == 0)
            {
                if(!noStore)
                    manager = ExpConfig::Setup::GetLastFound()->GetCalibrationDataManager();
                nCh = result.TaggEffs.size();
                vS.resize(nCh);
                vSy.resize(nCh);
                setupName = result.Setup;
            }
            if (n_TaggEffs > 0 && result.Setup != setupName )
                throw runtime_error("Different Setupnames within file list found!");

            for ( auto ch = 0u ; ch < nCh ; ++ch)
            {
                if (result.BkgFitChi2.at(ch) < chi2cut_channels )
                {
                    auto errorSqr = sqr(result.TaggEffErrors.at(ch));
                    vS.at(ch)   += 1.0 / errorSqr;
                    vSy.at(ch)  += result.TaggEffs.at(ch) / errorSqr;
                }
                else
                {
                    LOG(WARNING) << "Skipping channel " << ch << " for " << record.at(1) << ": bkg-chi2 = " << result.BkgFitChi2.at(ch) << ".";
                }
            }

            n_TaggEffs++;
        }
    }
    taggEff_t tagF(setupName,TID(),nCh);
    for ( auto ch = 0u; ch < nCh ; ++ch)
    {
        tagF.TaggEffs.at(ch)      = 1.0 * vSy.at(ch) / vS.at(ch);
        tagF.TaggEffErrors.at(ch) = sqrt(1.0 / vS.at(ch));
    }


    fillHistSingle(tagF.TaggEffs,tagF.TaggEffErrors);
    if (!noStore)
        storeResult(tagF,manager,calibration::TaggEff::GetDataName());
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
    if ( !graphTime )
    {
        graphTime = new TGraphErrors();
    }
    std_ext::RMS mte;
    std_ext::RMS mtee;

    for ( auto ch = 0u ; ch < data.size() ; ++ch )
    {
        hist_channels_2d->Fill(std_ext::to_iso8601(time).c_str(), ch,data.at(ch));
        mte.Add(data.at(ch));
        hist_channels_2d_errors->Fill(std_ext::to_iso8601(time).c_str(), ch,std::abs( 1.0 * dataErrors.at(ch) / data.at(ch)));
        mtee.Add(dataErrors.at(ch));
    }
    auto N = graphTime->GetN();
    graphTime->SetPoint(N,time,mte.GetMean());
    graphTime->SetPointError(N,0,mtee.GetMean());

}

void storeHist( const taggEff_t& result)
{
    if (histfac == nullptr)
        histfac = make_shared<HistogramFactory>(result.Setup);

    auto nCh = result.TaggEffs.size();

    auto hist = histfac->makeTH1D(result.Setup,"channel","#eta",BinSettings(nCh));
    for ( auto ch = 0u ; ch < nCh ; ++ch )
    {
        hist->Fill(ch,result.TaggEffs.at(ch));
        hist->SetBinError(ch+1,result.TaggEffErrors.at(ch));
    }
}

void storeResult(const taggEff_t& result, shared_ptr<calibration::DataManager> manager, const string& calibrationName,
                 const Calibration::AddMode_t& addMode)
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

    manager->Add(cdata, addMode);
}


