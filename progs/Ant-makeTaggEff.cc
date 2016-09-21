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
#include "base/WrapTFile.h"

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


struct channelHist_t
{

    TH1D* Hist = nullptr;
    shared_ptr<HistogramFactory> HistFac;
    string Title;

    channelHist_t(const string& title, shared_ptr<HistogramFactory> histfac):
        HistFac(histfac),
        Title(title) {}
    void Fill(const vector<double>& data, const vector<double>& dataErrors)
    {
        if ( !Hist )
            Hist = HistFac->makeTH1D(Title,"channel","Efficiency",BinSettings(data.size()));

        for ( auto ch = 0u ; ch < data.size() ; ++ch)
        {
            Hist->Fill(ch,data.at(ch));
            Hist->SetBinError(ch+1,dataErrors.at(ch)); // bin-numbering ( bin 0 is underflow ) ...
        }
    }
};

struct channelHistTime_t
{

    TH2D* Hist = nullptr;
    TH2D* HistErrors = nullptr;
    TGraphErrors* GraphMeans = nullptr;
    shared_ptr<HistogramFactory> HistFac;
    string Title;

    channelHistTime_t(const string& title, shared_ptr<HistogramFactory> histfac):
        HistFac(histfac),
        Title(title) {}
    void Fill(const vector<double>& data, const vector<double>& dataErrors, const unsigned time )
    {
        if ( !Hist )
            Hist = HistFac->makeTH2D(Title,
                                     "","channel",
                                     BinSettings(1,0,0),BinSettings(data.size()));
        if ( !HistErrors )
            HistErrors = HistFac->makeTH2D(std_ext::formatter() << Title << " relative errors",
                                           "","channel",
                                           BinSettings(1,0,0),BinSettings(data.size()));
        if ( !GraphMeans )
            GraphMeans = HistFac->make<TGraphErrors>();

        std_ext::RMS mte;
        std_ext::RMS mtee;
        for ( auto ch = 0u ; ch < data.size() ; ++ch )
        {
            Hist->Fill(std_ext::to_iso8601(time).c_str(), ch,data.at(ch));
            mte.Add(data.at(ch));
            HistErrors->Fill(std_ext::to_iso8601(time).c_str(), ch,std::abs( 1.0 * dataErrors.at(ch) / data.at(ch)));
            mtee.Add(dataErrors.at(ch));
        }
        auto N = GraphMeans->GetN();
        GraphMeans->SetPoint(N,time,mte.GetMean());
        GraphMeans->SetPointError(N,0,mtee.GetMean());
    }
};


static volatile bool interrupt = false;
static bool noStore = false;
static bool histOut = false;

auto failExit = [] (const string& message)
{
    LOG(ERROR) << message;
    exit(EXIT_FAILURE);
};

void storeClibrationData(const taggEff_t& result, shared_ptr<calibration::DataManager> manager, const string& calibrationName,
                         const Calibration::AddMode_t& addMode = Calibration::AddMode_t::RightOpen);
taggEffTriple_t* processFiles(const vector<string>& files, shared_ptr<channelHist_t> chHist);
taggEff_t mediateCSV(const vector<string>& csvFiles);
void processCSV(const string& csvFile, shared_ptr<channelHistTime_t> chHistTime);
taggEff_t processManualData(const string& dataFile, const string& setupName);

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
    auto cmd_output     = cmd.add<TCLAP::ValueArg<string>>("",  "save-hist", "Save results to a histogram", false, "","filename");

    //switches
    auto cmd_batchmode  = cmd.add<TCLAP::SwitchArg>("b", "batch",     "Run in batch mode (no ROOT shell afterwards)");
    auto cmd_nostore    = cmd.add<TCLAP::SwitchArg>("n", "nostore",   "don't store, only show results");

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
    histOut = cmd_output->isSet();

    unique_ptr<WrapTFileOutput> masterFile;
    if(histOut) {
        masterFile = std_ext::make_unique<WrapTFileOutput>(cmd_output->getValue(),
                                                    WrapTFileOutput::mode_t::recreate,
                                                     true);
    }

    auto histfac = make_shared<HistogramFactory>("makeTaggEff");

    shared_ptr<channelHistTime_t> chHistCsv = nullptr;
    if (mode_csv->isSet())
    {
        if (fileList.size() != 1)
            failExit("Provide one csv file!");
        chHistCsv = make_shared<channelHistTime_t>(fileList.at(0),histfac);
        processCSV(fileList.at(0),chHistCsv);
    }

    shared_ptr<channelHist_t> chHistMeanCsv;
    if (mode_csv_mean->isSet())
    {
        if (cmd_output)
        if (fileList.size() < 1)
            failExit("Provide at least one csv file!");
        auto result = mediateCSV(fileList);
        chHistMeanCsv = make_shared<channelHist_t>(result.Setup,histfac);
        chHistMeanCsv->Fill(result.TaggEffs,result.TaggEffErrors);

    }

    shared_ptr<channelHist_t> chHistGroup;
    taggEffTriple_t* triple_grp = nullptr;
    if (mode_group->isSet())
    {
        if (fileList.size() != 3)
            failExit("Group should contain three files in following order: 1st background, TaggEff-run, background 2.");
        chHistGroup = make_shared<channelHist_t>(fileList.at(1),histfac);
        triple_grp = processFiles(fileList,chHistGroup);
    }

    shared_ptr<channelHist_t> chHistManual;
    if (mode_manual->isSet())
    {
        if (!cmd_setup->isSet())
            failExit("Must provide setup name for data.");
        if (fileList.size() != 1)
            failExit("Provide one data file!");
        auto result = processManualData(fileList.at(0),cmd_setup->getValue());
        chHistManual = make_shared<channelHist_t>(std_ext::formatter() << result.Setup << " - manual data",histfac);
        chHistManual->Fill(result.TaggEffs,result.TaggEffErrors);
    }



    // OUPUT ==============================

    argc=1; // prevent TRint to parse any cmdline except prog name
    auto app = cmd_batchmode->isSet() || !std_ext::system::isInteractive()
               ? nullptr
               : std_ext::make_unique<TRint>("Ant-makeSigmas",&argc,argv,nullptr,0,true);
    if(app) {

        for ( auto i: {chHistMeanCsv,chHistGroup,chHistManual})
        {
            if (i)
            {
                i->Hist->SetStats(false);
                canvas("TaggEff") << drawoption("E") << i->Hist << endc;
            }
        }
        if (chHistCsv)
        {
            canvas("TaggEff") << drawoption("colz") << chHistCsv->Hist << endr
                              << chHistCsv->HistErrors << endr << endc;
            canvas("means") << drawoption("AP") << chHistCsv->GraphMeans << endc;
            chHistCsv->GraphMeans->GetXaxis()->SetTimeDisplay(1);
            chHistCsv->GraphMeans->GetXaxis()->SetTimeFormat("%d.%m.%Y %F 1970-01-01 00:00:00");
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

        if(masterFile)
            LOG(INFO) << "Stopped running, but close ROOT properly to write data to disk.";
        app->Run(kTRUE); // really important to return...
        if(masterFile)
            LOG(INFO) << "Writing output file...";
        masterFile = nullptr;   // and to destroy the master WrapTFile before TRint is destroyed
    }

    return EXIT_SUCCESS;
}

taggEff_t processManualData(const string& dataFile, const string& setupName)
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

    if (!noStore)
    {
        auto manager = ExpConfig::Setup::GetLastFound()->GetCalibrationDataManager();
        storeClibrationData(result,manager,calibration::TaggEff::GetDataName(),Calibration::AddMode_t::AsDefault);
    }

    return result;
}

taggEffTriple_t* processFiles(const vector<string>& files, shared_ptr<channelHist_t> chHist)
{
    auto taggEff = new taggEffTriple_t(files.at(0),files.at(1),files.at(2));
    taggEff_t result = taggEff->GetTaggEffSubtracted();
    auto manager = ExpConfig::Setup::GetLastFound()->GetCalibrationDataManager();

    chHist->Fill(result.TaggEffs,result.TaggEffErrors);

    if (!noStore)
        storeClibrationData(result,manager,calibration::TaggEff::GetDataName());

    return taggEff;
}

void processCSV(const string& csvFile, shared_ptr<channelHistTime_t> chHistTime)
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

        chHistTime->Fill(result.TaggEffs,result.TaggEffErrors,result.FirstID.Timestamp);
        if (!noStore)
            storeClibrationData(result,manager,calibration::TaggEff::GetDataName());
        setupName = result.Setup;
        n_TaggEffs++;
    }
}

taggEff_t mediateCSV(const vector<string>& csvFiles)
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
    if (!noStore)
        storeClibrationData(tagF,manager,calibration::TaggEff::GetDataName(),Calibration::AddMode_t::AsDefault);
    return tagF;
}

void storeClibrationData(const taggEff_t& result, shared_ptr<calibration::DataManager> manager, const string& calibrationName,
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


