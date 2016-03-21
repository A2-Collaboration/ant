#include "base/Logger.h"
#include "base/CmdLine.h"
#include <iostream>
#include "base/WrapTFile.h"
#include "TTree.h"
#include "tree/TID.h"

using namespace std;
using namespace ant;

//using cmd_t  = std::function<void(TTree*)>;
//using cmds_t = std::map<const std::string, cmd_t>;

TTree* getTIDTree(WrapTFileInput& file) {
    TTree* tree =nullptr;

    if(!file.GetObject("data_tid", tree))
        file.GetObject("h12_tid", tree);

    return tree;
}

int main (int argc, char** argv)
{

    SetupLogger();
    TCLAP::CmdLine cmd("Ant-calib - Fit histograms and calculate new calibration parameters", ' ', "0.1");
    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"level");
    auto cmd_treename = cmd.add<TCLAP::ValueArg<string>>("t","tree","Tree name", false, "","");
    auto cmd_cmd = cmd.add<TCLAP::ValueArg<string>>("c","cmd","command", true, "","");
    auto cmd_inputfiles  = cmd.add<TCLAP::UnlabeledMultiArg<string>>("inputfiles","Ant files with histograms",true,"inputfiles");
    cmd.parse(argc, argv);

    if(cmd_verbose->isSet())
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());

    if(cmd_cmd->getValue() == "tid_check") {

        if(cmd_inputfiles->getValue().size() != 2) {
            LOG(ERROR) << "exactly two files required for tid check!";
            return EXIT_FAILURE;
        }

        WrapTFileInput file0(cmd_inputfiles->getValue().at(0));
        WrapTFileInput file1(cmd_inputfiles->getValue().at(1));

        auto tree0 = getTIDTree(file0);

        if(!tree0) {
            LOG(ERROR) << "No TID tree in " << cmd_inputfiles->getValue().at(0);
            return 10;
        }

        auto tree1 = getTIDTree(file1);
        if(!tree1) {
            LOG(ERROR) << "No TID tree in " << cmd_inputfiles->getValue().at(1);
            return 11;
        }

        ant::TID* tid0 = nullptr;
        ant::TID* tid1 = nullptr;

        tree0->SetBranchAddress("tid", &tid0);
        tree1->SetBranchAddress("tid", &tid1);

        for(Long64_t i=0; i<10; ++i) {
            tree0->GetEntry(i);
            tree1->GetEntry(i);
            if(*tid0 != *tid1) {
                LOG(ERROR) << "TID missmatch!";
                return 20;
            }
        }

        LOG(INFO) << "All OK";

    } else {

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
    }

    return EXIT_SUCCESS;
}
