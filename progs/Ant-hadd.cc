#include "base/CmdLine.h"
#include "base/Logger.h"
#include "base/std_ext/string.h"
#include "base/std_ext/memory.h"

#include "TDirectory.h"
#include "TFile.h"
#include "TList.h"
#include "TKey.h"
#include "TClass.h"
#include "TH1.h"
#include "TFileMergeInfo.h"
#include "root-addons/analysis_codes/hstack.h"

#include <list>
#include <string>
#include <map>

using namespace std;
using namespace ant;

template<typename T>
using unique_ptrs_t = vector<unique_ptr<T>>;

using sources_t = unique_ptrs_t<const TDirectory>;


void MergeRecursive(TDirectory& target, const sources_t& sources)
{
    LOG(INFO) << target.GetPath();


    map<string, sources_t> dirs;

    map<string, unique_ptrs_t<TH1>>    hists;
    map<string, unique_ptrs_t<hstack>> stacks;

    for(auto& source : sources) {
        TList* keys = source->GetListOfKeys();
        if(!keys)
            continue;

        // first create a unique list of names,
        // this prevents object with different cycles
        TIter nextk(keys);
        string prev_keyname;
        while(auto key = dynamic_cast<TKey*>(nextk()))
        {
            const string keyname = key->GetName();
            if(prev_keyname == keyname)
                continue;
            prev_keyname = keyname;

            auto cl = TClass::GetClass(key->GetClassName());

            if(cl->InheritsFrom(TDirectory::Class())) {
                auto dir = dynamic_cast<TDirectory*>(key->ReadObj());
                dirs[keyname].emplace_back(dir);
            }
            else if(cl->InheritsFrom(TH1::Class())) {
                auto obj = dynamic_cast<TH1*>(key->ReadObj());
                hists[keyname].emplace_back(obj);
            }
            else if(cl->InheritsFrom(hstack::Class())) {
                auto obj = dynamic_cast<hstack*>(key->ReadObj());
                stacks[keyname].emplace_back(obj);
            }
        }
    }



    for(const auto& it_dirs : dirs) {
        auto newdir = target.mkdir(it_dirs.first.c_str());
        MergeRecursive(*newdir, it_dirs.second);
    }

    target.cd();
    TFileMergeInfo info(addressof(target)); // for calling Merge

    for(const auto& it_hists : hists) {
        auto& hists = it_hists.second;
        auto& first = hists.front();
        for(auto it = next(hists.begin()); it != hists.end(); ++it) {
            first->Add(it->get());
        }
        target.WriteTObject(first.get());
    }


    for(const auto& it_stacks : stacks) {
        auto& stacks = it_stacks.second;
        auto& first = stacks.front();
        TList c;
        for(auto it = next(stacks.begin()); it != stacks.end(); ++it) {
            c.Add(it->get());
        }
        first->Merge(addressof(c), addressof(info));
        target.WriteTObject(first.get());
    }
}

//___________________________________________________________________________
int main( int argc, char **argv )
{

   SetupLogger();

   TCLAP::CmdLine cmd("Ant-hadd - Merge ROOT objects in files", ' ', "0.1");
   auto cmd_filenames  = cmd.add<TCLAP::UnlabeledMultiArg<string>>("files","ROOT files, first one is output",true,"ROOT files");
   cmd.parse(argc, argv);

   list<string> filenames(cmd_filenames->getValue().begin(), cmd_filenames->getValue().end());

   for(auto filename :filenames) {
       if(std_ext::string_starts_with(filename, "-")) {
           LOG(ERROR) << "Found '" << filename << "' with starting - parsed as filename, might be wrongly spelled option. "
                      << "Prepend ./ to use it as inputfile.";
           return 1;
       }
   }

   if(filenames.size()<2) {
       LOG(ERROR) << "Provide at least an outputfile and one input file";
       exit(EXIT_FAILURE);
   }

   auto outputfile = std_ext::make_unique<TFile>(filenames.front().c_str(), "RECREATE");
   filenames.pop_front();

   sources_t sources;
   for(const auto& filename : filenames) {
       sources.emplace_back(std_ext::make_unique<TFile>(filename.c_str(), "READ"));
   }

   MergeRecursive(*outputfile, sources);

   LOG(INFO) << "Finished, writing file...";

   outputfile->Write();

   exit(EXIT_SUCCESS);
}
