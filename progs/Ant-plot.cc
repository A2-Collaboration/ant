#include "base/Logger.h"
#include "tclap/CmdLine.h"
#include "tclap/ValuesConstraintExtra.h"
#include "base/WrapTFile.h"
#include "base/std_ext/string.h"
#include "base/std_ext/system.h"
#include "base/ProgressCounter.h"
#include "analysis/physics/Plotter.h"

#include "tree/TAntHeader.h"
#include "expconfig/ExpConfig.h"

#include "TSystem.h"
#include "TRint.h"

#include <list>

using namespace ant;
using namespace ant::analysis;
using namespace std;


struct Plotter_list_entry {
    unique_ptr<Plotter> plotter;
    long long entries;
    Plotter_list_entry(unique_ptr<Plotter> p): plotter(std::move(p)), entries(plotter->GetNumEntries()) {}

    bool operator<(const Plotter_list_entry& other) const noexcept {
        return entries < other.entries;
    }
};
using plotter_list_t = std::list<Plotter_list_entry>;


volatile static bool interrupt = false;

int main(int argc, char** argv) {
    SetupLogger();

    signal(SIGINT, [] (int) { interrupt = true;} );



    TCLAP::CmdLine cmd("plot", ' ', "0.1");
    auto cmd_input  = cmd.add<TCLAP::ValueArg<string>>("i","input", "Input file", true,"","input");
    auto cmd_output = cmd.add<TCLAP::ValueArg<string>>("o","output","Output file",false,"","filename");

    TCLAP::ValuesConstraintExtra<decltype(ExpConfig::Setup::GetNames())> allowedsetupnames(ExpConfig::Setup::GetNames());
    auto cmd_setup  = cmd.add<TCLAP::ValueArg<string>>("s","setup","Choose setup manually by name",false,"", &allowedsetupnames);
    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"int");

    TCLAP::ValuesConstraintExtra<decltype(analysis::PlotterRegistry::GetList())> allowedPlotters(PlotterRegistry::GetList());
    auto cmd_plotters = cmd.add<TCLAP::MultiArg<string>>("p","Plotter","Plotter classes to run", true, &allowedPlotters);

    auto cmd_batchmode = cmd.add<TCLAP::MultiSwitchArg>("b","batch","Run in batch mode (no ROOT shell afterwards)",false);
    auto cmd_maxevents = cmd.add<TCLAP::ValueArg<int>>("m","maxevents","Process only max events",false,0,"maxevents");

    auto cmd_options = cmd.add<TCLAP::MultiArg<string>>("O","options","Options for all physics classes, key=value",false,"");

    cmd.parse(argc, argv);
    if(cmd_verbose->isSet()) {
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());
    }

    WrapTFileInput inputfile(cmd_input->getValue());

    // check if there's a previous AntHeader present,
    // which could tell us the SetupName
    if(!cmd_setup->isSet()) {
        TAntHeader* previous_AntHeader;
        if(inputfile.GetObject<TAntHeader>("AntHeader",previous_AntHeader)) {
            const auto& setupname = previous_AntHeader->SetupName;
            if(!setupname.empty()) {
                ExpConfig::Setup::SetByName(setupname);
                LOG(INFO) << "Setup name set to '" << setupname << "' from input file";
            }
            else
                LOG(WARNING) << "Found AntHeader in input files, but SetupName was empty";
        }
    }
    else {
        // override the setup name from cmd line
        ExpConfig::Setup::SetByName(cmd_setup->getValue());
        LOG(INFO) << "Commandline override setup name to '" << cmd_setup->getValue() << "'";
    }

    unique_ptr<WrapTFileOutput> masterFile;
    if(cmd_output->isSet()) {
        masterFile = std_ext::make_unique<WrapTFileOutput>(cmd_output->getValue(), true, WrapTFileOutput::mode_t::recreate);
    }

    plotter_list_t plotters;
    long long maxEntries = 0;
    {
        auto popts = make_shared<OptionsList>();

        if(cmd_options->isSet()) {
            for(const auto& opt : cmd_options->getValue()) {
                popts->SetOption(opt);
            }
        }

        for(const auto& plotter_name : cmd_plotters->getValue()) {
            try {
                plotters.emplace_back(PlotterRegistry::Create(plotter_name, inputfile, popts));
                maxEntries = max(maxEntries, plotters.back().entries);
            } catch(const exception& e) {
                LOG(ERROR) << "Could not create plotter \"" << plotter_name << "\": " << e.what();
                return EXIT_FAILURE;
            }
        }

        auto unused_popts = popts->GetUnused();
        if(!unused_popts.empty()) {
            LOG(ERROR) << "These plotter options where not recognized: " << unused_popts;
            LOG(INFO)  << "Did you mean: " << popts->GetNotFound();
            return EXIT_FAILURE;
        }
    }


    if(cmd_maxevents->isSet()) {
        maxEntries = min(maxEntries, static_cast<long long>(cmd_maxevents->getValue()));
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

    plotters.sort(); // sort by max entries

    auto p = plotters.begin();

    const auto advp = [&p,&plotters] (const long long& i) {
        while(i>=p->entries) {
            ++p;
            if(p==plotters.end())
                return false;
        }
        return true;
    };

    for(entry = 0; !interrupt && advp(entry) && entry < maxEntries; ++entry) {

        for(auto plotter = p; plotter!=plotters.end(); ++plotter) {
                plotter->plotter->ProcessEntry(entry);
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
