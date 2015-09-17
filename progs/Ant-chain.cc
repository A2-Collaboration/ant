/**
  * @file Ant-chain.cc
  * @brief hadd like tool to concatenate TTrees in ROOT files using TChains.
  *
  *        Creates a TChain in the output file for each TTree found in the first input file
  *        and uses AddFile() to add each input file to each TChain.
  */

//Ant
#include "base/Logger.h"
#include "base/CmdLine.h"
#include "base/WrapTFile.h"

//ROOT
#include "TDirectory.h"
#include "TChain.h"
#include "TTree.h"

//stl
#include <string>
#include <list>

using namespace std;
using namespace ant;

list<string> GetTreeNames(const std::vector<string> filenames) {

    list<string> chains;

    if(filenames.size() > 0) {
        WrapTFileInput testfile(filenames.at(0));

        const auto trees = testfile.GetListOf<TTree>();

        for(const auto tree : trees) {
            LOG(INFO) << "Found TTree " << tree->GetName();
            chains.emplace_back(string(tree->GetName()));
        }
    }

    return chains;
}

int main(int argc, char** argv) {
    SetupLogger();

    TCLAP::CmdLine cmd("Ant-chain", ' ', "0.1");
    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"level");
    auto cmd_output     = cmd.add<TCLAP::ValueArg<string>>("o","output","Output file",true,"","filename");
    auto cmd_inputfiles = cmd.add<TCLAP::UnlabeledMultiArg<string>>("inputfiles","input root files",true,"inputfiles");

    cmd.parse(argc, argv);


    if(cmd_verbose->isSet())
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());

    if(cmd_inputfiles->isSet()) {
        const auto& inputs = cmd_inputfiles->getValue();

        const auto chain_names = GetTreeNames(inputs);

        WrapTFileOutput outfile(cmd_output->getValue(),
                                WrapTFileOutput::mode_t::recreate,
                                true);
        list<TChain*> chains;
        for(const auto& name : chain_names) {
            auto chain = new TChain(name.c_str());
            chains.push_back(chain);
            gDirectory->Add(chain);
        }

        for(const auto& file : inputs) {
            for(auto chain : chains) {
                const auto res = chain->AddFile(file.c_str());
                if(res != 1) {
                    LOG(WARNING) << "Problem with " <<chain->GetName() << " and file " << file << " (" << res << ")";
                }
            }
            LOG(INFO) << "Added " << file;
        }

        LOG(INFO) << "Done";
        for(auto chain : chains) {
            LOG(INFO) << chain->GetName() << ": " << chain->GetEntries() << " Entries";
        }

    }

    return EXIT_SUCCESS;
}
