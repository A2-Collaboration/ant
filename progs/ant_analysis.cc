

#include "analysis/input/DataReader.h"
#include "analysis/input/ant/AntReader.h"
#include "analysis/input/goat/GoatReader.h"
#include "analysis/OutputManager.h"

#include "analysis/physics/Physics.h"
#include "analysis/physics/omega/omega.h"
#include "analysis/physics/common/DataOverview.h"
#include "analysis/physics/common/CandidatesAnalysis.h"

#include "expconfig/ExpConfig.h"

#include "base/std_ext.h"
#include "base/Logger.h"
#include "base/CmdLine.h"
#include "base/ReadTFiles.h"

#include "TRint.h"

#include <sstream>
#include <string>

using namespace std;
using namespace ant::output;
using namespace ant;
using namespace ant::analysis;

int main(int argc, char** argv) {
    SetupLogger();


    TCLAP::CmdLine cmd("Omega Analysis", ' ', "0.1");
    auto input_type = cmd.add<TCLAP::SwitchArg>("g","goat","Input is GoAT files",false);
    auto input  = cmd.add<TCLAP::MultiArg<string>>("i","input","ant root file (events)",true,"string");
    auto cmdline_setup  = cmd.add<TCLAP::ValueArg<string>>("","setup","Choose setup",false,"","string");
    auto cmdline_calibrations  = cmd.add<TCLAP::MultiArg<string>>("c","calibration","Calibrations to run",false,"string");
    auto physicsclasses  = cmd.add<TCLAP::MultiArg<string>>("p","physics","Physics Class to run",false,"string");
    auto output = cmd.add<TCLAP::ValueArg<string>>("o","output","Output file",false,"","string");
    auto batchmode = cmd.add<TCLAP::SwitchArg>("b","batch","Run in batch mode (No ROOT Windows)",false);
    auto verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"int");

    cmd.parse(argc, argv);
    if(verbose->isSet()) {
        el::Loggers::setVerboseLevel(verbose->getValue());
    }

    ExpConfig::ManualSetupName = cmdline_setup->getValue(); // force setupname
    using setup_ptr_t = shared_ptr<ExpConfig::Setup>;
    setup_ptr_t setup = ExpConfig::Setup::Get(ExpConfig::ManualSetupName);
    if(setup != nullptr) {
        LOG(INFO) << "Found setup " << setup->GetName();
    }

    list<shared_ptr<Calibration::PhysicsModule>> enabled_calibrations;
    if(cmdline_calibrations->isSet()) {
        if(setup==nullptr) {
            stringstream ss_setups;
            for(const setup_ptr_t& item : ExpConfig::Setup::GetAll()) {
                ss_setups << item->GetName() << " ";
            }
            LOG(INFO)  << "Available setups: " << ss_setups.str();
            LOG(ERROR) << "Please specify a --setup if you want to use calibrations as physics modules";
        }
        else {
            stringstream ss_calibrations;

            for(const auto& calibration : setup->GetCalibrations()) {
                ss_calibrations << calibration->GetName() << " ";
                if(!std_ext::contains(cmdline_calibrations->getValue(), calibration->GetName()))
                    continue;
                enabled_calibrations.emplace_back(move(calibration));
            }
            if(enabled_calibrations.empty()) {
                LOG(ERROR) << "No calibrations enabled. Available: " << ss_calibrations.str();
            }
        }
    }


    PhysicsManager pm;
    OutputManager om;

    if(output->isSet())
        om.SetNewOutput(output->getValue());

    for(const auto& classname : physicsclasses->getValue()) {
        try {
            pm.AddPhysics( PhysicsRegistry::Create(classname) );
            VLOG(3) << "Activated physics class \"" << classname << "\"";
        } catch (...) {
            LOG(WARNING) << "Physics class \"" << classname << "\" is not registered";
        }
    }

    for(const auto& calibration : enabled_calibrations) {
        pm.AddPhysics(calibration->GetPhysicsModule());
    }

    // build the general ROOT file manager first
    auto filemanager = make_shared<ReadTFiles>();
    for(const auto& inputfile : input->getValue()) {
        VLOG(5) << "ROOT File Manager: Looking at file " << inputfile;
        if(filemanager->OpenFile(inputfile))
            LOG(INFO) << "Opened file '" << inputfile << "' as ROOT file";
        else
            VLOG(5) << "Could not add " << inputfile << " to ROOT file manager";
    }

    std::unique_ptr<input::DataReader> reader;

    if(input_type->getValue()) {
        reader = std_ext::make_unique<input::GoatReader>(filemanager);
    } else {
        reader = std_ext::make_unique<input::AntReader>(filemanager);
    }

    pm.ReadFrom(*reader.get());

    if(!batchmode->isSet()) {
        int a=0;
        char** b=nullptr;
        TRint app("omega",&a,b);
        pm.ShowResults();
        app.Run(kTRUE);
    }

    return 0;

}

