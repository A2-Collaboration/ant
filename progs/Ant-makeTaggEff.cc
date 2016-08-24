
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

string processFiles(const string& bkg1, const string& run, const string& bkg2)
{
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
            ant::TAntHeader* h;
            wrapFile.GetObject("AntHeader",h);
            setupName = h->SetupName;
        }
    };

    treeContainer_t tree_bkg1(bkg1);
    treeContainer_t tree_bkg2(bkg2);
    treeContainer_t tree_run(run);

    if ( !(tree_bkg1.setupName == tree_bkg2.setupName &&
         tree_bkg2.setupName == tree_run.setupName ) )
    {
        LOG(ERROR) << "Files in TaggEff-triple not from same Setup!";
        exit(EXIT_FAILURE);
    }

    return tree_bkg1.setupName; // any is ok it is checked
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
