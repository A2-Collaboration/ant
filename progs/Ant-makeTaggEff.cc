#include <string>
#include <map>

#include "analysis/physics/common/ProcessTaggEff.h"
#include "analysis/plot/RootDraw.h"

#include "tclap/CmdLine.h"
#include "base/Detector_t.h"
#include "base/Logger.h"
#include "base/PlotExt.h"
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
            GraphMeans = HistFac->makeGraphErrors("means","means");

        std_ext::RMS mte;
        std_ext::RMS mtee;
        for ( auto ch = 0u ; ch < data.size() ; ++ch )
        {
            Hist->Fill(std_ext::to_iso8601(time).c_str(), ch,data.at(ch));
            mte.Add(data.at(ch));
            HistErrors->Fill(std_ext::to_iso8601(time).c_str(), ch,std::abs( 1.0 * dataErrors.at(ch) / data.at(ch)));
            mtee.Add(dataErrors.at(ch));
        }
        GraphExt::FillGraphErrors(GraphMeans,time,mte.GetMean(),0.,mtee.GetMean());
    }
};


// chi2 cut value used for background fit to determine if a channel should be skipped
static double chi2cut_channels;
// some other bools used in various methods
static volatile bool interrupt = false;
static bool noStore = false;
static bool histOut = false;

auto failExit = [] (const string& message)
{
    LOG(ERROR) << message;
    exit(EXIT_FAILURE);
};

void storeCalibrationData(const taggEff_t& result, Calibration::AddMode_t addMode);
taggEffTriple_t* processFiles(const vector<string>& files, shared_ptr<channelHist_t> chHist,const HistogramFactory& histfac);
taggEff_t mediateCSV(const vector<string>& csvFiles,const HistogramFactory& histfac);
void processCSV(const string& csvFile, shared_ptr<channelHistTime_t> chHistTime,const HistogramFactory& histfac);

int main( int argc, char** argv )
{
    SetupLogger();

    signal(SIGINT, [] (int) {
        LOG(INFO) << ">>> Interrupted";
        interrupt = true;
    });

    TCLAP::CmdLine cmd("Ant-makeTaggEff - Create Tagging Efficiencies from triples background1 - run - background2", ' ', "0.1");

    // implemented modes
    TCLAP::SwitchArg mode_csv("", "csv", "CSV file with comma-separated groups of bkg1,run,bkg2 - calibration data will be added as right open patches for each group");
    TCLAP::SwitchArg mode_csv_mean("", "csv-mean", "CSV file with comma-separated groups of bkg1,run,bkg2 - calibration data will be the mean over all measurements");
    TCLAP::SwitchArg mode_group("", "group", "Provide a single group containing 'background1 run background2' files for adding a right open patch");
    // register modes here:
    vector<TCLAP::Arg*> modes;
    modes.push_back(&mode_csv);
    modes.push_back(&mode_csv_mean);
    modes.push_back(&mode_group);
    // add mutually exclusive modes as argument
    cmd.xorAdd(modes);

    // other settings:
    auto cmd_output    = cmd.add<TCLAP::ValueArg<string>>("o", "output", "Output file", false, "", "filename");
    auto cmd_chi2      = cmd.add<TCLAP::ValueArg<double>>("c", "chi2", "chi2 value used as cut condition to skip channels while fitting the background", false, 15.0, "chi2 cut value");

    // switches
    auto cmd_batchmode = cmd.add<TCLAP::SwitchArg>("b", "batch",   "Run in batch mode (no ROOT shell afterwards)");
    auto cmd_nostore   = cmd.add<TCLAP::SwitchArg>("n", "nostore", "Don't store Tagging Efficiencies in the calibration database, only show results");

    // files
    auto cmd_filelist  = cmd.add<TCLAP::UnlabeledMultiArg<string>>("inputfiles", "Input files to read from", true, "inputfiles");


    cmd.parse(argc, argv);

    auto fileList = cmd_filelist->getValue();
    noStore = cmd_nostore->isSet();
    histOut = cmd_output->isSet();
    chi2cut_channels = cmd_chi2->getValue();

    unique_ptr<WrapTFileOutput> masterFile;
    if (histOut) {
        masterFile = std_ext::make_unique<WrapTFileOutput>(cmd_output->getValue(), true);
    }

    auto histfac = make_shared<HistogramFactory>("makeTaggEff");

    shared_ptr<channelHistTime_t> chHistCsv = nullptr;
    if (mode_csv.isSet())
    {
        if (fileList.size() != 1)
            failExit("Provide one CSV file!");
        chHistCsv = make_shared<channelHistTime_t>(fileList.at(0),histfac);
        processCSV(fileList.at(0),chHistCsv,*histfac);
    }

    shared_ptr<channelHist_t> chHistMeanCsv;
    if (mode_csv_mean.isSet())
    {
        if (cmd_output)
        if (fileList.size() < 1)
            failExit("Provide at least one CSV file!");
        auto result = mediateCSV(fileList,*histfac);
        chHistMeanCsv = make_shared<channelHist_t>(result.Setup,histfac);
        chHistMeanCsv->Fill(result.TaggEffs,result.TaggEffErrors);

    }

    shared_ptr<channelHist_t> chHistGroup;
    taggEffTriple_t* triple_grp = nullptr;
    if (mode_group.isSet())
    {
        if (fileList.size() != 3)
            failExit("Group should contain three files in following order: 1st background, TaggEff-run, 2nd background");
        chHistGroup = make_shared<channelHist_t>(fileList.at(1),histfac);
        triple_grp = processFiles(fileList,chHistGroup,*histfac);
    }


    // OUTPUT ==============================

    argc=1; // prevent TRint to parse any cmdline except prog name
    auto app = cmd_batchmode->isSet() || !std_ext::system::isInteractive()
               ? nullptr
               : std_ext::make_unique<TRint>("Ant-makeSigmas",&argc,argv,nullptr,0,true);
    if(app) {

        for ( auto i: {chHistMeanCsv,chHistGroup})
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
            canvas control("control");
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

taggEffTriple_t* processFiles(const vector<string>& files, shared_ptr<channelHist_t> chHist, const HistogramFactory& histfac)
{
    auto taggEff = new taggEffTriple_t(files.at(0),files.at(1),files.at(2),histfac);
    taggEff_t result = taggEff->GetTaggEffSubtracted();

    chHist->Fill(result.TaggEffs,result.TaggEffErrors);

    if (!noStore) {
        storeCalibrationData(result, Calibration::AddMode_t::RightOpen);
    }

    return taggEff;
}

void processCSV(const string& csvFile, shared_ptr<channelHistTime_t> chHistTime, const HistogramFactory& histfac)
{
    auto histLambda     = histfac.makeTH1D("Decay constants","decay constant [1/s]","#",BinSettings(100,0,0),"histDecayConst");
    auto graphLambda    = histfac.makeGraph("Decay Constants","graphDecayConsts");

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

        taggEffTriple_t taggEff(record.at(0),record.at(1),record.at(2),histfac);
        taggEff_t result = taggEff.GetTaggEffSubtracted();
        histLambda->Fill(taggEff.GetDecayConstant());
        GraphExt::FillGraph(graphLambda,result.FirstID.Timestamp,taggEff.GetDecayConstant());


        // in case of first tagging efficiency (which doesn't have to be at the first day of the beamtime)
        // take the first date of this specific beamtime specified via its Setup
        if (n_TaggEffs == 0) {
            ExpConfig::Setup::SetByName(result.Setup);
            auto start = ExpConfig::Setup::Get().GetStartDate();
            result.FirstID = TID(to_time_t(start, false));
        }

        if (n_TaggEffs > 0 && result.Setup != setupName )
                throw runtime_error("Different Setupnames within file list found!");

        chHistTime->Fill(result.TaggEffs,result.TaggEffErrors,result.FirstID.Timestamp);
        if (!noStore) {
            storeCalibrationData(result, Calibration::AddMode_t::RightOpen);
        }
        setupName = result.Setup;
        n_TaggEffs++;
    }
}

taggEff_t mediateCSV(const vector<string>& csvFiles, const HistogramFactory& histfac)
{
    string setupName;
    size_t n_TaggEffs(0);
    size_t nCh(0);

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

            const auto result = taggEffTriple_t(record.at(0),record.at(1),record.at(2),histfac).GetTaggEffSubtracted();

            if (n_TaggEffs == 0)
            {
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
    if (!noStore) {
        storeCalibrationData(tagF, Calibration::AddMode_t::AsDefault);
    }
    return tagF;
}

void storeCalibrationData(const taggEff_t& result, Calibration::AddMode_t addMode)
{
    auto tagger = ExpConfig::Setup::GetDetector<TaggerDetector_t>();
    auto manager = ExpConfig::Setup::Get().GetCalibrationDataManager();

    const auto& calibrationName = calibration::TaggEff::GetModuleName(tagger->Type);

    TCalibrationData cdata(
                calibrationName,
                result.FirstID,
                result.FirstID // its right open interval (or default)
                );

    for ( auto ch=0u ; ch < result.TaggEffs.size() ; ++ch )
    {
        cdata.Data.emplace_back(ch,result.TaggEffs.at(ch));
        // see calibration::TaggEff module: errors stored in fit-parameters!!!
        cdata.FitParameters.emplace_back(ch,vector<double>(1,result.TaggEffErrors.at(ch)));
    }

    manager->Add(cdata, addMode);
}


