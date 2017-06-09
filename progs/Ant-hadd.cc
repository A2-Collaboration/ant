#include "tclap/CmdLine.h"
#include "base/Logger.h"
#include "base/std_ext/string.h"
#include "base/std_ext/memory.h"
#include "base/ProgressCounter.h"
#include "base/std_ext/system.h"

#include "root-addons/analysis_codes/hadd.h"

#include "TSystem.h"
#include "TROOT.h"
#include "TInterpreter.h"
#include "TFileMerger.h"
#include "TFile.h"
#include "TClass.h"

#include <list>
#include <string>
#include <map>

using namespace std;
using namespace ant;

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
   hadd::sources_t sources;
   for(const auto& filename : filenames) {
       sources.emplace_back(std_ext::make_unique<TFile>(filename.c_str(), "READ"));
   }

   // progress updates only when running interactively
   if(std_ext::system::isInteractive())
       ProgressCounter::Interval = 3;

   unsigned nPaths = 0;
   ProgressCounter progress([&nPaths] (chrono::duration<double> elapsed) {
       LOG(INFO) << nPaths/elapsed.count() << " paths/s";
       nPaths = 0;
   });

   hadd::MergeRecursive(*outputfile, sources, nPaths);

   LOG(INFO) << "Finished, writing file " << outputfile->GetName();

   outputfile->Write();

   exit(EXIT_SUCCESS);
}
