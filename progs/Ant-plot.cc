#include "base/Logger.h"

#include "analysis/plot/CutTree.h"
#include "analysis/physics/etaprime/etaprime_omega_gamma.h"

#include "expconfig/ExpConfig.h"
#include "expconfig/detectors/EPT.h"
#include <list>
#include "base/CmdLine.h"
#include "base/interval.h"
#include "base/printable.h"
#include "base/WrapTFile.h"
#include "base/std_ext/string.h"
#include "base/std_ext/system.h"
#include "base/ProgressCounter.h"
#include "analysis/physics/Plotter_Traits.h"

#include "TSystem.h"
#include "TRint.h"

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::plot;
using namespace std;

volatile bool interrupt = false;


int main(int argc, char** argv) {
    SetupLogger();

    signal(SIGINT, [] (int) { interrupt = true; cerr << "Stopping..." << endl; } );

    TCLAP::CmdLine cmd("plot", ' ', "0.1");
    auto cmd_input  = cmd.add<TCLAP::ValueArg<string>>("i","input", "Input file", true,"","input");
    auto cmd_output = cmd.add<TCLAP::ValueArg<string>>("o","output","Output file",false,"","filename");
    auto cmd_plotters = cmd.add<TCLAP::MultiArg<string>>("p","Plotter","Plotter classes to run",true,"plotter");

    auto cmd_batchmode = cmd.add<TCLAP::MultiSwitchArg>("b","batch","Run in batch mode (no ROOT shell afterwards)",false);
    auto cmd_maxevents = cmd.add<TCLAP::ValueArg<int>>("m","maxevents","Process only max events",false,0,"maxevents");


    auto cmd_setupname = cmd.add<TCLAP::ValueArg<string>>("s","setup","Override setup name", false, "Setup_2014_07_EPT_Prod", "setup");

    cmd.parse(argc, argv);

    const auto setup_name = cmd_setupname->getValue() ;
    auto setup = ExpConfig::Setup::Get(setup_name);
    if(setup == nullptr) {
        LOG(ERROR) << "Did not find setup instance for name " << setup_name;
        return 1;
    }

    WrapTFileInput input(cmd_input->getValue());

    unique_ptr<WrapTFileOutput> masterFile;
    if(cmd_output->isSet()) {
        masterFile = std_ext::make_unique<WrapTFileOutput>(cmd_output->getValue(),
                                                    WrapTFileOutput::mode_t::recreate,
                                                     true); // cd into masterFile upon creation
    }

    list<pair<unique_ptr<Plotter_Trait>, long long>> plotters;

    OptionsPtr PlotterOpts = nullptr;

    long long maxEntries = cmd_maxevents->isSet() ? cmd_maxevents->getValue() : 0;

    for(const auto& plotter_name : cmd_plotters->getValue()) {
        try {
            plotters.emplace_back(PlotterRegistry::Create(plotter_name, input, PlotterOpts), -1);
            plotters.back().second = plotters.back().first->GetNumEntries();
            maxEntries = max(maxEntries, plotters.back().second);
        } catch(exception& e) {
            LOG(INFO) << "Could not create plotter \"" << plotter_name << "\": " << e.what();
        }
    }

    for(long long entry = 0; entry < maxEntries && !interrupt; ++entry) {

        for(auto& plotter : plotters) {
            if(entry < plotter.second) {    //@todo: make this more efficient
                plotter.first->ProcessEntry(entry);
            }
        }

        if(interrupt)
            break;
    }

    for(auto& plotter : plotters) {
        plotter.first->Finish();
    }


    if(!cmd_batchmode->isSet()) {
        if(!std_ext::system::isInteractive()) {
            LOG(INFO) << "No TTY attached. Not starting ROOT shell.";
        }
        else {

            argc=0; // prevent TRint to parse any cmdline
            TRint app("EtapOmegaG_plot",&argc,argv,nullptr,0,true);

            if(masterFile)
                LOG(INFO) << "Stopped running, but close ROOT properly to write data to disk.";

            for(auto& plotter : plotters) {
                plotter.first->ShowResult();
            }

            app.Run(kTRUE); // really important to return...
            if(masterFile)
                LOG(INFO) << "Writing output file...";
            masterFile = nullptr;   // and to destroy the master WrapTFile before TRint is destroyed
        }
    }

    return EXIT_SUCCESS;
}
