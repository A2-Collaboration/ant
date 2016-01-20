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
#include "base/std_ext/string.h"
#include "base/std_ext/system.h"

//ROOT
#include "TDirectory.h"
#include "TChain.h"
#include "TTree.h"

//stl
#include <string>
#include <list>
#include <functional>

using namespace std;
using namespace ant;

template<typename T>
string get_path(T* dir) {
    if(dir && dir->GetMotherDir()) {
        auto prev = get_path(dir->GetMotherDir());
        return std_ext::formatter() << prev << (prev.empty() ? "" : "/") << dir->GetName();
    }
    else {
        return "";
    }
}

set<string> GetTreeNames(const std::vector<string> filenames) {
    set<string> chains; // we use a set to avoid adding the TTree more than once
    if(filenames.size() > 0) {
        WrapTFileInput firstfile(filenames.front());
        firstfile.Traverse([&chains] (TKey* key) {
            if(string(key->GetClassName()) == "TTree")
                chains.insert(get_path(key));
        });
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
        if(chain_names.empty()) {
            LOG(ERROR) << "No TTree found in first input file";
            return EXIT_FAILURE;
        }

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
                const auto absFile = std_ext::system::absolutePath(file);
                const auto res = chain->AddFile(absFile.c_str());
                if(res != 1) {
                    LOG(WARNING) << "Problem with " <<chain->GetName() << " and file " << file << " (" << res << ")";
                }
            }
            LOG(INFO) << "Added " << file;
        }

        for(auto chain : chains) {
            LOG(INFO) << chain->GetName() << ": " << chain->GetEntries() << " Entries";
            if(std_ext::contains(chain->GetName(), "/"))
                LOG(INFO) << "Remember to use 'TTree* tree = gDirectory->GetKey(\"" << chain->GetName() << "\")->ReadObj()'";
        }

    }

    return EXIT_SUCCESS;
}
