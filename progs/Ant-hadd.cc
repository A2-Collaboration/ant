#include "tclap/CmdLine.h"
#include "base/Logger.h"
#include "base/std_ext/string.h"
#include "base/std_ext/memory.h"
#include "base/ProgressCounter.h"
#include "base/std_ext/system.h"

#include "tree/TAntHeader.h"

#include "TDirectory.h"
#include "TFile.h"
#include "TList.h"
#include "TKey.h"
#include "TClass.h"
#include "TH1.h"
#include "TFileMergeInfo.h"
#include "root-addons/analysis_codes/hstack.h"

#include "TSystem.h"
#include "TROOT.h"
#include "TInterpreter.h"
#include "TFileMerger.h"


#include <list>
#include <string>
#include <map>

using namespace std;
using namespace ant;

template<typename T>
using unique_ptrs_t = vector<unique_ptr<T>>;

using sources_t = unique_ptrs_t<const TDirectory>;

unsigned nPaths = 0;

template<typename T>
struct pair_t {
    explicit pair_t(const string& name) : Name(name) {}
    string Name;
    T Item;
};

template<typename C, typename... Args>
void add_by_name(C& c, const string& name, Args&&... args) {
    using T = typename C::value_type;
    auto it = std::find_if(c.begin(), c.end(), [name] (const T& item) {
        return item.Name == name;
    });
    if(it == c.end()) {
        c.emplace_back(name);
        c.back().Item.emplace_back(std::forward<Args>(args)...);
    }
    else {
        it->Item.emplace_back(std::forward<Args>(args)...);
    }
}

void MergeRecursive(TDirectory& target, const sources_t& sources)
{
    nPaths++;
    ProgressCounter::Tick();

    vector<pair_t<sources_t>> dirs;

    vector<pair_t<unique_ptrs_t<TH1>>>    hists;
    vector<pair_t<unique_ptrs_t<hstack>>> stacks;
    vector<pair_t<unique_ptrs_t<TAntHeader>>> headers;

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
                add_by_name(dirs, keyname, dir);
            }
            else if(cl->InheritsFrom(TH1::Class())) {
                auto obj = dynamic_cast<TH1*>(key->ReadObj());
                add_by_name(hists, keyname, obj);
            }
            else if(cl->InheritsFrom(hstack::Class())) {
                auto obj = dynamic_cast<hstack*>(key->ReadObj());
                add_by_name(stacks, keyname, obj);
            }
            else if(cl->InheritsFrom(TAntHeader::Class())) {
                auto obj = dynamic_cast<TAntHeader*>(key->ReadObj());
                add_by_name(headers, keyname, obj);
            }
        }
    }



    for(const auto& it_dirs : dirs) {
        auto newdir = target.mkdir(it_dirs.Name.c_str());
        MergeRecursive(*newdir, it_dirs.Item);
    }

    target.cd();
    TFileMergeInfo info(addressof(target)); // for calling Merge

    for(const auto& it_hists : hists) {
        auto& items = it_hists.Item;

        // check if at least one hist has labels,
        // the others could be never filled (so ROOT treats them as normal hists)
        const auto hasLabels = [] (const unique_ptr<TH1>& h) {
            return h->GetXaxis()->GetLabels() != nullptr;
        };
        const auto it_h_withLabels = std::find_if(items.begin(), items.end(), hasLabels);

        auto& first = items.front();
        if(it_h_withLabels != items.end()) {

            // again, scan the histograms for empty hists without labels
            // IMHO, this is a bug in ROOT that empty hists cannot be merged with labeled hists
            auto& h_withLabels = *it_h_withLabels;
            for(auto& h : items) {
                if(hasLabels(h))
                    continue;
                for(int bin=0;bin<h->GetNbinsX()+1;bin++)
                    if(h->GetBinContent(bin) != 0)
                        throw std::runtime_error("Found non-empty unlabeled hist "
                                                 + string(h->GetDirectory()->GetPath()));
                // prepare the axis labels of the empty hist, labeled hist should have at least
                // one bin filled
                h->Fill(h_withLabels->GetXaxis()->GetBinLabel(1), 0.0);
            }

            auto& first = items.front();
            TList c;
            for(auto it = next(items.begin()); it != items.end(); ++it) {
                c.Add(it->get());
            }
            first->Merge(addressof(c));
        }
        else {
            for(auto it = next(items.begin()); it != items.end(); ++it) {
                first->Add(it->get());
            }
        }
        target.WriteTObject(first.get());
    }

    for(const auto& it : stacks) {
        auto& items = it.Item;
        auto& first = items.front();
        TList c;
        for(auto it = next(items.begin()); it != items.end(); ++it) {
            c.Add(it->get());
        }
        first->Merge(addressof(c), addressof(info));
        target.WriteTObject(first.get());
    }

    for(const auto& it : headers) {
        auto& items = it.Item;
        auto& first = items.front();
        TList c;
        for(auto it = next(items.begin()); it != items.end(); ++it) {
            c.Add(it->get());
        }
        first->Merge(addressof(c));
        target.WriteTObject(first.get());
    }
}

void do_nativemode(const string& outputfile, const list<string>& inputfiles) {
    Bool_t force = kTRUE; // changed from defaults
    Bool_t skip_errors = kFALSE;
    Bool_t reoptimize = kFALSE;
    Bool_t noTrees = kFALSE;
    Int_t maxopenedfiles = 0;
    Int_t verbosity = 99;
    Int_t newcomp = 1; // compression level


    gSystem->Load("libTreePlayer");
    TClass::GetClass("ROOT::Cintex::Cintex"); // autoload Cintex if it exist.
    if (gInterpreter->IsLoaded("libCintex")) {
       gROOT->ProcessLine("ROOT::Cintex::Cintex::Enable();");
    }

    TFileMerger merger(kFALSE,kFALSE);
    merger.SetMsgPrefix("Ant-hadd");
    merger.SetPrintLevel(verbosity - 1);
    if (maxopenedfiles > 0) {
       merger.SetMaxOpenedFiles(maxopenedfiles);
    }
    if (!merger.OutputFile(outputfile.c_str(),force,newcomp) ) {
       LOG(ERROR) << "hadd error opening target file";
       exit(EXIT_FAILURE);
    }

    for ( auto& inputfile : inputfiles ) {
        if( ! merger.AddFile(inputfile.c_str()) ) {
          if ( skip_errors ) {
             LOG(ERROR) << "hadd skipping file with error: " << inputfile << endl;
          } else {
             LOG(ERROR) << "hadd exiting due to error in " << inputfile << endl;
             exit(EXIT_FAILURE);
          }
       }
    }
    if (reoptimize) {
       merger.SetFastMethod(kFALSE);
    } else {
       if (merger.HasCompressionChange()) {
          // Don't warn if the user any request re-optimization.
          LOG(WARNING) <<"hadd Sources and Target have different compression levels"<<endl;
          LOG(WARNING) <<"hadd merging will be slower"<<endl;
       }
    }

    merger.SetNotrees(noTrees);
    Bool_t status = merger.Merge();

    if(!status) {
        LOG(ERROR) << "hadd failure during the merge of " << merger.GetMergeList()->GetEntries() << " input files in " << outputfile << ".\n";
        exit(EXIT_FAILURE);
    }

    LOG(INFO) << "hadd merged " << merger.GetMergeList()->GetEntries() << " input files in " << outputfile << ".\n";
}

//___________________________________________________________________________
int main( int argc, char **argv )
{

   SetupLogger();

   TCLAP::CmdLine cmd("Ant-hadd - Merge ROOT objects in files", ' ', "0.1");
   auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"int");
   auto cmd_nativemode = cmd.add<TCLAP::MultiSwitchArg>("","native","Run native TFileMerger, is slow on large trees",false);
   auto cmd_filenames  = cmd.add<TCLAP::UnlabeledMultiArg<string>>("files","ROOT files, first one is output",true,"ROOT files");
   cmd.parse(argc, argv);
   if(cmd_verbose->isSet()) {
       el::Loggers::setVerboseLevel(cmd_verbose->getValue());
   }

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

   // important to copy the std::string here
   const auto outputfilename = filenames.front();
   filenames.pop_front();

   if(cmd_nativemode->isSet()) {
       do_nativemode(outputfilename, filenames);
       exit(EXIT_SUCCESS);
   }

   auto outputfile = std_ext::make_unique<TFile>(outputfilename.c_str(), "RECREATE");
   sources_t sources;
   for(const auto& filename : filenames) {
       sources.emplace_back(std_ext::make_unique<TFile>(filename.c_str(), "READ"));
   }

   // progress updates only when running interactively
   if(std_ext::system::isInteractive())
       ProgressCounter::Interval = 3;

   ProgressCounter progress([] (chrono::duration<double> elapsed) {
       LOG(INFO) << nPaths/elapsed.count() << " paths/s";
       nPaths = 0;
   });

   MergeRecursive(*outputfile, sources);

   LOG(INFO) << "Finished, writing file " << outputfile->GetName();

   outputfile->Write();

   exit(EXIT_SUCCESS);
}
