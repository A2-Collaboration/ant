#include <map>

#include "analysis/plot/HistogramFactories.h"

#include "base/std_ext/system.h"
#include "base/Logger.h"
#include "base/CmdLine.h"
#include "base/WrapTFile.h"
#include "base/std_ext/string.h"

#include "expconfig/setups/Setup.h"
#include "expconfig/ExpConfig.h"

#include "analysis/plot/root_draw.h"

#include "TTree.h"
#include "TRint.h"
#include "TH1D.h"

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
            vector<double> electrons;
            double livetime = 0;
            vector<double> tdcs;
        };

        means_t getMeans() const
        {
            means_t means;
            means.electrons.resize(nchannels);
            means.tdcs.resize(nchannels);

            auto nentries = wrapTree.Tree->GetEntries();

            for (auto entry = 0 ; entry < nentries ; ++entry)
            {
                wrapTree.Tree->GetEntry(entry);
                for ( auto channel = 0u ; channel < wrapTree.TaggRates().size(); ++channel)
                {
                    means.electrons.at(channel) += 1.0 * wrapTree.TaggRates().at(channel);
                    means.tdcs.at(channel)      += 1.0 * wrapTree.TDCHits().at(channel);
                }
                means.livetime += 1.0 * wrapTree.ExpLivetime;

            }

            means.livetime = means.livetime / nentries;
            for (auto& e: means.electrons)
                e = e / nentries;
            for (auto& t: means.tdcs)
                t = t / nentries;

            return means;
        }

    };

    treeContainer_t bkg1;
    treeContainer_t run;
    treeContainer_t bkg2;



public:

    taggEffTriple_t(const string& bkg1f_, const string& runf_, const string& bkg2f_):
        bkg1(bkg1f_),
        run(runf_),
        bkg2(bkg2f_)
    {
        if ( !(bkg1.setupName == bkg2.setupName &&
             bkg2.setupName == run.setupName ) )
            throw runtime_error("Files in TaggEff-triple not from same Setup!");
    }

    string SetupName() const{return bkg1.setupName;}

    // startID for these taggEff values should be the end of this taggEff-triple:
    TID startID() const
    {
        bkg2.wrapTree.Tree->GetEntry(bkg2.wrapTree.Tree->GetEntries()-1);
        return bkg2.wrapTree.LastID();
    }

    // TODO: calibration data
    resultSet_t  GetTaggEff() const
    {
        resultSet_t result(SetupName(),startID());

        treeContainer_t::means_t m_bkg1 = bkg1.getMeans();
        treeContainer_t::means_t m_run = run.getMeans();
        treeContainer_t::means_t m_bkg2 = bkg2.getMeans();

        auto nchannels = bkg1.nchannels;
        result.TaggEffs.resize(nchannels,0);
        result.TaggEffErrors.resize(nchannels,0);

        for (auto channel = 0u ; channel < nchannels ; ++channel)
        {
            result.TaggEffs.at(channel)  =  (m_run.tdcs.at(channel)
                                            - ( (m_bkg1.tdcs.at(channel) + m_bkg2.tdcs.at(channel)) / 2.0 ) ) ;
            result.TaggEffs.at(channel) *= 1.0 / ( m_run.livetime * m_run.electrons.at(channel)
                                            - ( (m_bkg1.livetime * m_bkg1.electrons.at(channel)
                                                + m_bkg2.livetime * m_bkg2.electrons.at(channel) ) / 2.0) );
        }


        return result;
    }

};


void processFiles(const string& bkg1, const string& run, const string& bkg2)
{
    taggEffTriple_t taggEff(bkg1,run,bkg2);
    resultSet_t result = taggEff.GetTaggEff();
    auto manager = ExpConfig::Setup::GetLastFound()->GetCalibrationDataManager();

    storeResult(result.FirstID,manager,calibration::TaggEff::GetDataName(),result.TaggEffs);
}

bool processCSV(const string& csvFile)
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
            storeResult(it->second,manager,calibration::TaggEff::GetDataName(),result.TaggEffs);
        }
        if (n_TaggEffs > 0)
        {
            if(result.Setup != setupName )
                throw runtime_error("Different Setupnames within file list found!");
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
        if (processCSV(inCsv))
            return EXIT_SUCCESS;
    }

    const auto inBkg1    = cmd_bkg1->getValue();
    const auto inTaggEff = cmd_taggEff->getValue();
    const auto inBkg2    = cmd_bkg2->getValue();

    processFiles(inBkg1,inTaggEff,inBkg2);

    return EXIT_SUCCESS;
}
