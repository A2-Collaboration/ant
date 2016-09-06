#include <map>

#include "analysis/plot/HistogramFactories.h"

#include "base/std_ext/system.h"
#include "base/Logger.h"
#include "base/CmdLine.h"
#include "base/WrapTFile.h"
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

#include "calibration/DataManager.h"
#include "tree/TCalibrationData.h"
#include "calibration/modules/TaggEff.h"

#include "analysis/physics/common/ProcessTaggEff.h"
#include "base/Detector_t.h"
#include "expconfig/detectors/EPT.h"

#include "tree/TAntHeader.h"


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

struct resultSet_t
{
    string          Setup;
    TID             FirstID;
    vector<double>  TaggEffs;
    vector<double>  TaggEffErrors;
    resultSet_t(const string& setup, const TID& firstID):
        Setup(setup),
        FirstID(firstID),
        TaggEffs(),
        TaggEffErrors() {}
};

void fillHistSingle( const vector<double>& data)
{
    if ( !hist_channels )
    {
        hist_channels = new TH1D("taggEff", "Generated Tagging Efficiencies", data.size(),0,data.size());
        hist_channels->SetXTitle("channel");
        hist_channels->SetYTitle("Efficiency");
    }
    for ( auto ch = 0u ; ch < data.size() ; ++ch)
        hist_channels->Fill(ch,data.at(ch));
}

void fillHistTime( const vector<double>& data, const unsigned time )
{
    if ( !hist_channels_2d )
    {
        hist_channels_2d = new TH2D("taggEffTime","Generated Tagging Efficiencies", 1 ,0,0,data.size(),0,data.size());
        hist_channels_2d->SetXTitle("");
        hist_channels_2d->SetYTitle("channel");
    }

    for ( auto ch = 0u ; ch < data.size() ; ++ch )
        hist_channels_2d->Fill(std_ext::to_iso8601(time).c_str(), ch,data.at(ch));
}

void storeResult(const TID& startID, shared_ptr<calibration::DataManager> manager, const string& calibrationName, vector<double> data)
{
    TCalibrationData cdata(
                calibrationName,
                startID,
                startID
                );

    for ( auto ch=0u ; ch < data.size() ; ++ch )
        cdata.Data.emplace_back(ch,data.at(ch));

    manager->Add(cdata, Calibration::AddMode_t::RightOpen);
}

class taggEffTriple_t
{
protected:
    struct treeContainer_t
    {
        WrapTFileInput wrapFile;
        physics::ProcessTaggEff::TreeScalarReads wrapTree;

        string setupName;
        size_t nchannels;


        treeContainer_t(const string& filename):
            wrapFile(filename)
        {
            auto treeName = physics::ProcessTaggEff::treeAccessName();
            if (!wrapFile.GetObject(treeName,wrapTree.Tree))
            {
                LOG(ERROR) << "Cannot find tree " << treeName
                           << " in file " << filename;
                exit(EXIT_FAILURE);
            }

            wrapTree.LinkBranches(wrapTree.Tree);

            ant::TAntHeader* h;
            wrapFile.GetObject("AntHeader",h);
            setupName = h->SetupName;
            ExpConfig::Setup::SetManualName(setupName);
            nchannels = ExpConfig::Setup::GetDetector<TaggerDetector_t>()->GetNChannels();
        }

        struct means_t {
            vector<double> Scalers;
            vector<double> ScalersErrors;
            double Livetime = 0;
            double LivetimeError = 0;
            vector<double> Tdcs;
            vector<double> TdcsErrors;
        };

        means_t getMeans() const
        {
            means_t means;
            means.Scalers.resize(nchannels);
            means.ScalersErrors.resize(nchannels);
            means.Tdcs.resize(nchannels);
            means.TdcsErrors.resize(nchannels);

            vector<std_ext::RMS> rms_scalers(nchannels);
            vector<std_ext::RMS> rms_tds(nchannels);
            std_ext::RMS rms_livetime;


            auto nentries = wrapTree.Tree->GetEntries();

            for (auto entry = 0 ; entry < nentries ; ++entry)
            {
                wrapTree.Tree->GetEntry(entry);
                for ( auto channel = 0u ; channel < nchannels; ++channel)
                {
                    rms_scalers.at(channel).Add(1.0 * wrapTree.TaggRates().at(channel));
                    rms_tds.at(channel).Add(1.0 * wrapTree.TDCRates().at(channel));
                }
                rms_livetime.Add(1.0 * wrapTree.ExpLivetime);
            }

            for ( auto channel = 0u ; channel < nchannels ; ++channel)
            {
                means.Scalers.at(channel)       = rms_scalers.at(channel).GetMean();
                means.ScalersErrors.at(channel) = rms_scalers.at(channel).GetRMS();

                means.Tdcs.at(channel)          = rms_tds.at(channel).GetMean();
                means.TdcsErrors.at(channel)    = rms_tds.at(channel).GetRMS();
            }
            means.Livetime      = rms_livetime.GetMean();
            means.LivetimeError = rms_livetime.GetRMS();

            return means;
        }

    };

    treeContainer_t bkg1;
    treeContainer_t run;
    treeContainer_t bkg2;

    TID startID;


public:

    taggEffTriple_t(const string& bkg1f_, const string& runf_, const string& bkg2f_):
        bkg1(bkg1f_),
        run(runf_),
        bkg2(bkg2f_)
    {
        LOG(INFO) << "Loading TaggEff measurement for:  "
                  << "1st bkg: " << bkg1f_ << "  "
                  << "run: "     << runf_  << "  "
                  << "2nd bkg: " << bkg2f_;

        if ( !(bkg1.setupName == bkg2.setupName &&
             bkg2.setupName == run.setupName ) )
            throw runtime_error("Files in TaggEff-triple not from same Setup!");
        WrapTFileInput file(bkg2f_);
        ant::TAntHeader* header;
        file.GetObject("AntHeader",header);
        if (!header)
            throw runtime_error("No Ant header in found!");
        startID = header->LastID;
    }

    string SetupName() const{return bkg1.setupName;}

   unsigned sanityChecks(const treeContainer_t::means_t& bkg1,
                 const treeContainer_t::means_t& run,
                 const treeContainer_t::means_t& bkg2) const
    {
        auto severity = 0u;
        auto nchannels = bkg1.Scalers.size();




        // per file checks:
        for ( auto m: {bkg1, run, bkg2})
        {
            for (auto ch = 0u ; ch < nchannels ; ++ch)
            {
                if (m.Scalers.at(ch) < m.Tdcs.at(ch))
                {
                    LOG(ERROR) << "Detected file with higher TDC rate than scaler rate!";
                    LOG(ERROR) << m.Scalers.at(ch) << " < " << m.Tdcs.at(ch);
                    severity++;
                }
            }
            if ( !(0 <= m.Livetime && m.Livetime <= 1))
            {
                LOG(ERROR) << "Detected file with live time not in [0,1]!";
                LOG(ERROR) << "livetime: " << m.Livetime;
                severity++;
            }
        }

        // bkg
        for ( auto ch = 0u ; ch < nchannels ; ++ch)
        {
            if ( run.Scalers.at(ch) < bkg1.Scalers.at(ch) || run.Scalers.at(ch) < bkg2.Scalers.at(ch))
                LOG(WARNING) << "Detected file with higher background than actual measurement in scalers, channel " << ch << "!";
            if ( run.Tdcs.at(ch) < bkg1.Tdcs.at(ch) || run.Tdcs.at(ch) < bkg2.Tdcs.at(ch))
                LOG(WARNING) << "Detected file with higher background than actual measurement in TDCs, channel " << ch << "!";
        }
        return severity;
    }



    const resultSet_t  GetTaggEff() const
    {
        resultSet_t result(SetupName(),startID);

        auto denum = [] (const double L,  const double s,
                          const double L1, const double s1,
                          const double L2, const double s2)
        {
            return L * s - ( ( L1 * s1 + L2 * s2 ) / 2.0 );
        };
        auto denum_bkg = [] (const double L,  const double s,
                          const double L1, const double s1,
                          const double L2, const double s2)
        {
            return 2 * L *s - L1 * s1 - L2 * s2;
        };


        const auto m_bkg1 = bkg1.getMeans();
        const auto m_run = run.getMeans();
        const auto m_bkg2 = bkg2.getMeans();
        if ( sanityChecks(m_bkg1,m_run,m_bkg2) > 0 )
            throw runtime_error("Tagg-Eff-triple didn't pass sanity checks, doublecheck input files!");

        auto nchannels = bkg1.nchannels;
        result.TaggEffs.resize(nchannels,0);
        result.TaggEffErrors.resize(nchannels,0);

        for (auto ch = 0u ; ch < nchannels ; ++ch)
        {
            auto denum_ch         = denum(m_run.Livetime,  m_run.Scalers.at(ch),
                                          m_bkg1.Livetime, m_bkg1.Scalers.at(ch),
                                          m_bkg2.Livetime, m_bkg2.Scalers.at(ch));
            auto denum_bkg_ch_sqr = denum_bkg(m_run.Livetime,  m_run.Scalers.at(ch),
                                          m_bkg1.Livetime, m_bkg1.Scalers.at(ch),
                                          m_bkg2.Livetime, m_bkg2.Scalers.at(ch));


            result.TaggEffs.at(ch)  = 1.0 * m_run.Tdcs.at(ch) / denum_ch;

            result.TaggEffErrors.at(ch) =
                    sqrt(   sqr( m_run.TdcsErrors.at(ch)/denum_ch )
                          + sqr( m_run.LivetimeError * m_run.Tdcs.at(ch) * m_run.Scalers.at(ch)  / sqr(denum_ch) )
                          + sqr( m_run.ScalersErrors.at(ch) * m_run.Livetime * m_run.Tdcs.at(ch) / sqr(denum_ch) )
                          + sqr( m_bkg1.LivetimeError * 2.0 * m_run.Tdcs.at(ch) * m_bkg1.Scalers.at(ch)  / denum_bkg_ch_sqr )
                          + sqr( m_bkg1.ScalersErrors.at(ch) * 2.0 * m_run.Tdcs.at(ch) * m_bkg1.Livetime / denum_bkg_ch_sqr )
                          + sqr( m_bkg2.LivetimeError * 2.0 * m_run.Tdcs.at(ch) * m_bkg2.Scalers.at(ch)  / denum_bkg_ch_sqr )
                          + sqr( m_bkg2.ScalersErrors.at(ch) * 2.0 * m_run.Tdcs.at(ch) * m_bkg2.Livetime / denum_bkg_ch_sqr )
                         );

            if ( !(IntervalD(0,1).Contains(result.TaggEffs.at(ch))) )
                LOG(WARNING)  << "Calculated Tagging efficiency for channel " << ch
                              << " not in [0,1]: "
                              << result.TaggEffs.at(ch) << "!";
        }

        return result;
    }

};


void processFiles(const string& bkg1, const string& run, const string& bkg2, const bool noStore)
{
    taggEffTriple_t taggEff(bkg1,run,bkg2);
    resultSet_t result = taggEff.GetTaggEff();
    auto manager = ExpConfig::Setup::GetLastFound()->GetCalibrationDataManager();

    fillHistSingle(result.TaggEffs);

    if (!noStore)
        storeResult(result.FirstID,manager,calibration::TaggEff::GetDataName(),result.TaggEffs);
}

bool processCSV(const string& csvFile, const bool noStore)
{
    shared_ptr<calibration::DataManager> manager = nullptr;

    ifstream csvStream(csvFile);
    if (!csvStream)
    {
        LOG(ERROR) << "Error reading File list " << csvFile << ".";
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
        resultSet_t result = taggEff.GetTaggEff();

        //check if setup is valid for this method --> String in map?
        auto it = startIDs.find(result.Setup);
        if (it == startIDs.end())
            throw runtime_error("Setup not valid for csv mode!");

        if (n_TaggEffs == 0)
        {
            manager = ExpConfig::Setup::GetLastFound()->GetCalibrationDataManager();
            fillHistTime(result.TaggEffs,it->second.Timestamp);
            if (!noStore)
                storeResult(it->second,manager,calibration::TaggEff::GetDataName(),result.TaggEffs);
        }
        if (n_TaggEffs > 0)
        {
            if(result.Setup != setupName )
                throw runtime_error("Different Setupnames within file list found!");
            fillHistTime(result.TaggEffs,result.FirstID.Timestamp);
            if (!noStore)
                storeResult(result.FirstID,manager,calibration::TaggEff::GetDataName(),result.TaggEffs);
        }
        setupName = result.Setup;
        n_TaggEffs++;
    }

    return true;
}


int main( int argc, char** argv )
{
    SetupLogger();

    signal(SIGINT, [] (int) {
        LOG(INFO) << ">>> Interrupted";
        interrupt = true;
    });

    TCLAP::CmdLine cmd("Ant-makeTaggEff", ' ', "0.1");

    auto cmd_input        = cmd.add<TCLAP::ValueArg<string>>("i","input",      "CSV file with bkg1-run-bkg2 groups",             false,  "", "csv");
    auto cmd_bkg1         = cmd.add<TCLAP::ValueArg<string>>("f","first-bkg",  "input for first background",                     false,  "", "root");
    auto cmd_taggEff      = cmd.add<TCLAP::ValueArg<string>>("t","taggEff",    "TaggEff - run",                                  false,  "", "root");
    auto cmd_bkg2         = cmd.add<TCLAP::ValueArg<string>>("s","second-bkg", "input for second background",                    false,  "", "root");

    auto cmd_batchmode    = cmd.add<TCLAP::MultiSwitchArg>  ("b","batch",      "Run in batch mode (no ROOT shell afterwards)",   false);
    auto cmd_nostore      = cmd.add<TCLAP::SwitchArg>("","nostore","don't store, only show results");

    cmd.parse(argc, argv);

    const auto inCsv   = cmd_input->getValue();

    if (cmd_input->isSet())
    {
        if ( cmd_bkg1->isSet()    ||
             cmd_bkg2->isSet()    ||
             cmd_taggEff->isSet()    )
        {
            LOG(ERROR) << "Provide either csv file list or triple of backgrounds and TaggEff run!";
            return EXIT_FAILURE;
        }
        processCSV(inCsv, cmd_nostore->isSet());
    }

    if (!cmd_input->isSet())
    {
        const auto inBkg1    = cmd_bkg1->getValue();
        const auto inTaggEff = cmd_taggEff->getValue();
        const auto inBkg2    = cmd_bkg2->getValue();

        processFiles(inBkg1,inTaggEff,inBkg2,cmd_nostore->isSet());
    }
    argc=1; // prevent TRint to parse any cmdline except prog name
    auto app = cmd_batchmode->isSet() || !std_ext::system::isInteractive()
               ? nullptr
               : std_ext::make_unique<TRint>("Ant-makeSigmas",&argc,argv,nullptr,0,true);

    if(app) {
        canvas c("TaggEff");
        if (hist_channels)
            c << hist_channels;
        if (hist_channels_2d)
            c << drawoption("colz") << hist_channels_2d;
        c << endc;

        app->Run(kTRUE); // really important to return...
    }


    return EXIT_SUCCESS;
}
