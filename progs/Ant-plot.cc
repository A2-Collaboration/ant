#include "base/Logger.h"
#include "base/CmdLine.h"
#include "base/detail/tclap/ValuesConstraintExtra.h"
#include "base/printable.h"
#include "base/WrapTFile.h"
#include "base/std_ext/string.h"
#include "base/std_ext/system.h"
#include "base/ProgressCounter.h"
#include "analysis/physics/Plotter.h"

#include "TSystem.h"
#include "TRint.h"

#include <list>

using namespace ant;
using namespace ant::analysis;
using namespace std;


struct Plotter_list_entry {
    unique_ptr<Plotter> plotter;
    long long entries;
    Plotter_list_entry(unique_ptr<Plotter> p): plotter(std::move(p)), entries(p->GetNumEntries()) {}
};
using plotter_list_t = std::list<Plotter_list_entry>;


volatile static bool interrupt = false;

int main(int argc, char** argv) {
    SetupLogger();

    signal(SIGINT, [] (int) { interrupt = true;} );



    TCLAP::CmdLine cmd("plot", ' ', "0.1");
    auto cmd_input  = cmd.add<TCLAP::ValueArg<string>>("i","input", "Input file", true,"","input");
    auto cmd_output = cmd.add<TCLAP::ValueArg<string>>("o","output","Output file",false,"","filename");

    TCLAP::ValuesConstraintExtra<decltype(analysis::PlotterRegistry::GetList())> allowedPlotters(PlotterRegistry::GetList());
    auto cmd_plotters = cmd.add<TCLAP::MultiArg<string>>("p","Plotter","Plotter classes to run", true, &allowedPlotters);

    auto cmd_batchmode = cmd.add<TCLAP::MultiSwitchArg>("b","batch","Run in batch mode (no ROOT shell afterwards)",false);
    auto cmd_maxevents = cmd.add<TCLAP::ValueArg<int>>("m","maxevents","Process only max events",false,0,"maxevents");

    cmd.parse(argc, argv);



    WrapTFileInput input(cmd_input->getValue());

    unique_ptr<WrapTFileOutput> masterFile;
    if(cmd_output->isSet()) {
        masterFile = std_ext::make_unique<WrapTFileOutput>(cmd_output->getValue(), true, WrapTFileOutput::mode_t::recreate);
    }

    plotter_list_t plotters;

    OptionsPtr PlotterOpts = nullptr;

    long long maxEntries = cmd_maxevents->isSet() ? cmd_maxevents->getValue() : 0;

    for(const auto& plotter_name : cmd_plotters->getValue()) {
        try {
            plotters.emplace_back(PlotterRegistry::Create(plotter_name, input, PlotterOpts));
        } catch(exception& e) {
            LOG(INFO) << "Could not create plotter \"" << plotter_name << "\": " << e.what();
        }
    }

    if(plotters.empty()) {
        LOG(ERROR) << "No active plotters. Nothing to do.";
    }


    long long entry;

    ProgressCounter progress(
                [&entry, maxEntries]
                (std::chrono::duration<double> elapsed)
    {
        const double percent = double(entry)/maxEntries;

        static double last_PercentDone = 0;
        const double speed = (percent - last_PercentDone)/elapsed.count();
        LOG(INFO) << setw(2) << std::setprecision(4)
                  << percent*100 << " % done, ETA: " << ProgressCounter::TimeToStr((1-percent)/speed);
        last_PercentDone = percent;
    });

    ProgressCounter::Interval = 5; //sec

    for(entry = 0; entry < maxEntries && !interrupt; ++entry) {

        for(auto& plotter : plotters) {
            if(entry < plotter.entries) {    //@todo: make this more efficient
                plotter.plotter->ProcessEntry(entry);
            }
        }

        ProgressCounter::Tick();

        if(interrupt)
            break;
    }


    LOG(INFO) << "Analyzed " << entry << " records"
              << ", speed " << entry/progress.GetTotalSecs() << " event/s";

    for(auto& plotter : plotters) {
        plotter.plotter->Finish();
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
                plotter.plotter->ShowResult();
            }

            app.Run(kTRUE); // really important to return...
            if(masterFile)
                LOG(INFO) << "Writing output file...";

            //cleanup before TRint is destroyed:
            masterFile = nullptr;
            plotters.clear();
        }
    }

    return EXIT_SUCCESS;
}
