#include "base/Logger.h"
#include "base/CmdLine.h"
#include <iostream>
#include "base/WrapTFile.h"
#include "TTree.h"

using namespace std;
using namespace ant;

//using cmd_t  = std::function<void(TTree*)>;
//using cmds_t = std::map<const std::string, cmd_t>;


int main (int argc, char** argv)
{

    SetupLogger();
    TCLAP::CmdLine cmd("Ant-calib - Fit histograms and calculate new calibration parameters", ' ', "0.1");
    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"level");
    auto cmd_treename = cmd.add<TCLAP::ValueArg<string>>("t","tree","Tree name", true, "","");
    auto cmd_cmd = cmd.add<TCLAP::ValueArg<string>>("c","cmd","command", true, "","");
    auto cmd_inputfiles  = cmd.add<TCLAP::UnlabeledMultiArg<string>>("inputfiles","Ant files with histograms",true,"inputfiles");
    cmd.parse(argc, argv);

    if(cmd_verbose->isSet())
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());

    const auto treename = cmd_treename->getValue();

    for(const auto& file : cmd_inputfiles->getValue()) {

        WrapTFileInput infile;

        try {
            infile.OpenFile(file);
        } catch (std::exception& e) {
            LOG(WARNING) << e.what();
        }

        TTree* tree = nullptr;
        infile.GetObject(treename, tree);

        if(tree) {
            if(cmd_cmd->getValue() == "GetEntries") {
                cout << file << ": " << tree->GetEntries() << endl;
            } else {
                LOG(ERROR) << "unknown command";
            }
        }
    }

    return EXIT_SUCCESS;
}
