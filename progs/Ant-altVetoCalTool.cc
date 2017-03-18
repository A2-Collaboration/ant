/**
  * @file OmegaEtaG2_plot.cc
  * @brief Plotting tool for TTrees from the ant::analysis::physics::OmegaEtaG2 physics class
  */

#include "base/Logger.h"

#include "analysis/physics/pi0/JustPi0.h"
#include "analysis/utils/ParticleTools.h"

#include "tclap/CmdLine.h"
#include "base/interval.h"
#include "base/WrapTFile.h"
#include "base/std_ext/string.h"
#include "base/std_ext/system.h"
#include "base/std_ext/math.h"
#include "tree/TAntHeader.h"

#include <list>
#include <vector>
#include <algorithm>
#include <memory>

#include "TSystem.h"
#include "TRint.h"

#include "TStyle.h"
#include "TCutG.h"

using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace std;

volatile static bool interrupt = false;


int main(int argc, char** argv) {
    SetupLogger();

    signal(SIGINT, [] (int) { interrupt = true; } );

    TCLAP::CmdLine cmd("plot", ' ', "0.1");
    auto cmd_input = cmd.add<TCLAP::ValueArg<string>>("i","input","Input file",true,"","input");
    auto cmd_batchmode = cmd.add<TCLAP::MultiSwitchArg>("b","batch","Run in batch mode (no ROOT shell afterwards)",false);
    auto cmd_maxevents = cmd.add<TCLAP::MultiArg<int>>("m","maxevents","Process only max events",false,"maxevents");
    auto cmd_output = cmd.add<TCLAP::ValueArg<string>>("o","output","Output file",false,"","filename");

    auto cmd_tree = cmd.add<TCLAP::ValueArg<string>>("","tree","Tree name",false,"T1","treename");


    cmd.parse(argc, argv);


    WrapTFileInput input(cmd_input->getValue());

    auto link_branches = [&input] (const string treename, WrapTTree* wraptree, long long expected_entries) {
        TTree* t;
        if(!input.GetObject(treename,t))
            throw runtime_error("Cannot find tree "+treename+" in input file");
        if(expected_entries>=0 && t->GetEntries() != expected_entries)
            throw runtime_error("Tree "+treename+" does not have entries=="+to_string(expected_entries));
        if(wraptree->Matches(t,false)) {
            wraptree->LinkBranches(t);
            return true;
        }
        return false;
    };

    TAntHeader* header = nullptr;
    input.GetObject("AntHeader", header);

    if(!header)
        LOG(ERROR) << "No AntHeder found in input!";


    physics::JustPi0::MultiPi0::MultiPi0Tree tree;

    if(!link_branches("JustPi0/m2Pi0/tree", addressof(tree), -1)) {
        LOG(WARNING) << "Cannot link branches of tree";
       // return 1;
    }

    const auto entries = tree.Tree->GetEntries();

    unique_ptr<WrapTFileOutput> masterFile;
    if(cmd_output->isSet()) {
        // cd into masterFile upon creation
        masterFile = std_ext::make_unique<WrapTFileOutput>(cmd_output->getValue(), true);
    }

    masterFile->WriteObject(header, "AntHeader");


    HistogramFactory PID_HistFac("PID_Energy");
    HistogramFactory VETO_HistFac("TAPSVeto_Energy");

    auto h_pid_bananas = PID_HistFac.makeTH3D(
                "PID Bananas",
                "CB Energy / MeV",
                "PID Energy / MeV",
                "Channel",
                BinSettings(400,0,800),
                BinSettings(200,0,18),
                24,
                "Bananas"
                );

    auto h_veto_bananas = VETO_HistFac.makeTH3D(
                "TAPSVeto Bananas",
                "TAPS LG Energy / MeV",
                "TAPSVeto Energy / MeV",
                "Channel",
                BinSettings(400,0,800),
                BinSettings(200,0,18),
                384,
                "Bananas"
                );


    LOG(INFO) << "Tree entries=" << entries;

    auto max_entries = entries;
    if(cmd_maxevents->isSet() && cmd_maxevents->getValue().back()<entries) {
        max_entries = cmd_maxevents->getValue().back();
        LOG(INFO) << "Running until " << max_entries;
    }

    for(long long entry=0;entry<max_entries;entry++) {
        if(interrupt)
            break;

        tree.Tree->GetEntry(entry);

        // Fill hists
        TH3* h = nullptr;

        if(tree.proton_det == 1) {
            h = h_pid_bananas;
        } else if(tree.proton_det == 2) {
            h = h_veto_bananas;
        }

        if(h)
            h->Fill(tree.proton_fitted().E()-ParticleTypeDatabase::Proton.Mass(),
                    tree.proton_vetoE,
                    tree.proton_vetoCh,
                    tree.Tagg_W
                    );


        if(entry % 100000 == 0)
            LOG(INFO) << "Processed " << 100.0*entry/entries << " %";
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

            app.Run(kTRUE); // really important to return...
            if(masterFile)
                LOG(INFO) << "Writing output file...";
            masterFile = nullptr;   // and to destroy the master WrapTFile before TRint is destroyed
        }
    }


    return 0;
}
