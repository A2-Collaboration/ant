#include <assert.h>

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

#include "calibration/Editor.h"
#include "calibration/DataManager.h"

#include "analysis/physics/common/ProcessTaggEff.h"

#include "tree/TAntHeader.h"

using namespace ant;
using namespace std;
using namespace ant::analysis;
using namespace ant::std_ext;
static volatile bool interrupt = false;


class taggEffTriple_t
{
protected:
    struct treeContainer_t
    {
        WrapTFileInput wrapFile;
        physics::ProcessTaggEff::TreeScalarReads wrapTree;

        string setupName;

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
        }
    };

    treeContainer_t bkg1tree;
    treeContainer_t runtree;
    treeContainer_t bkg2tree;

public:

    taggEffTriple_t(const string& bkg1, const string& run, const string& bkg2):
        bkg1tree(bkg1),
        runtree(run),
        bkg2tree(bkg2)
    {

        if ( !(bkg1tree.setupName == bkg2tree.setupName &&
             bkg2tree.setupName == runtree.setupName ) )
        {
            LOG(ERROR) << "Files in TaggEff-triple not from same Setup!";
            exit(EXIT_FAILURE);
        }

        ExpConfig::Setup::SetManualName(bkg1tree.setupName);

    }

    string const SetupName() const{ return bkg1tree.setupName;}

    // TODO: calibration data
    vector<double> const Get()
    {
        auto tagger = ExpConfig::Setup::GetDetector<TaggerDetector_t>();

        //found setup should match
        assert(tagger->GetNchannels() == bkg1tree.wrapTree.TDCHits().size());

//        for (int i = 0 ; i < tagger->GetNChannels() ; ++i )
//            calcTaggeff(i)

        return {};
    }

};


string processFiles(const string& bkg1, const string& run, const string& bkg2)
{

    taggEffTriple_t taggEff(bkg1,run,bkg2);

    cout << "TaggEff by channel:" << endl;
    for (const auto eff: taggEff.Get())
        cout << eff << "  ";
    cout << endl;




    return taggEff.SetupName(); // any is ok, it is checked
}

bool processCSV(const string& csvFile)
{
    ifstream csvStream(csvFile);
    if (!csvStream)
    {
        LOG(ERROR) << "Error reading File list " << csvFile << ".";
        return false;
    }

    string setupName;
    size_t n_TaggEffs(0);

    while (csvStream)
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
        {
            LOG(ERROR) << "Found line with wrong number of files, check your file list.";
            return false;
        }

        string currentSetup(processFiles(record[0],record[1],record[2]));
        if (n_TaggEffs > 0)
            if (currentSetup != setupName)
            {
                LOG(ERROR) << "Different Setupnames within file list found!";
                return false;
            }
        setupName = currentSetup;
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

    TCLAP::CmdLine cmd("Ant-makeSigmas", ' ', "0.1");

    auto cmd_input        = cmd.add<TCLAP::ValueArg<string>>("i","input",      "CSV file with bkg1-run-bkg2 groups",             false,  "", "csv");
    auto cmd_bkg1         = cmd.add<TCLAP::ValueArg<string>>("f","first-bkg",  "input for first background",                     false,  "", "root");
    auto cmd_taggEff      = cmd.add<TCLAP::ValueArg<string>>("t","taggEff",    "TaggEff - run",                                  false,  "", "root");
    auto cmd_bkg2         = cmd.add<TCLAP::ValueArg<string>>("s","second-bkg", "input for second background",                    false,  "", "root");


    auto cmd_output       = cmd.add<TCLAP::ValueArg<string>>("o","",           "output-file",                                    false, "", "rootfile");

    auto cmd_batchmode    = cmd.add<TCLAP::MultiSwitchArg>  ("b","batch",      "Run in batch mode (no ROOT shell afterwards)",   false);


    cmd.parse(argc, argv);

    const auto inCsv   = cmd_input->getValue();


    const auto outfile = cmd_output->getValue();


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
    const auto inBkg2    = cmd_bkg2->getValue();
    const auto inTaggEff = cmd_taggEff->getValue();

    LOG(INFO) << "Processed files for Setup " << processFiles(inBkg1,inBkg2,inTaggEff) << ".";

    return EXIT_SUCCESS;
}
