/**
  * @file Ant.cc
  * @brief Main Ant executable
  *
  */

#include "analysis/input/DataReader.h"
#include "analysis/input/ant/AntReader.h"
#include "analysis/input/goat/GoatReader.h"
#include "analysis/input/pluto/PlutoReader.h"
#include "analysis/utils/ParticleID.h"
#include "analysis/physics/PhysicsManager.h"

#include "expconfig/ExpConfig.h"
#include "expconfig/setups/SetupRegistry.h"

#include "calibration/DataManager.h"
#include "calibration/DataBase.h"

#include "unpacker/Unpacker.h"
#include "unpacker/RawFileReader.h"

#include "reconstruct/Reconstruct.h"

#include "tree/TAntHeader.h"

#include "base/WrapTFile.h"
#include "base/Logger.h"
#include "tclap/CmdLine.h"
#include "tclap/ValuesConstraintExtra.h"
#include "base/WrapTFile.h"
#include "base/std_ext/system.h"
#include "base/std_ext/container.h"
#include "base/GitInfo.h"

#include "TRint.h"
#include "TSystem.h"
#include "TROOT.h"

#include <sstream>
#include <string>
#include <csignal>

using namespace std;
using namespace ant;

volatile bool interrupt = false;
volatile bool terminated = false;


int main(int argc, char** argv) {
    SetupLogger();

    signal(SIGINT, [] (int) {
        cout << ">>> Interrupted" << endl;
        interrupt = true;
    });

    signal(SIGTERM, [] (int) {
        cout << ">>> Terminated" << endl;
        interrupt = true;
        terminated = true;
    });

    TCLAP::CmdLine cmd("Ant", ' ', "0.1");

    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"int");
    auto cmd_input  = cmd.add<TCLAP::MultiArg<string>>("i","input","Input files",true,"filename");

    TCLAP::ValuesConstraintExtra<decltype(ExpConfig::Setup::GetNames())> allowedsetupnames(ExpConfig::Setup::GetNames());
    auto cmd_setup  = cmd.add<TCLAP::ValueArg<string>>("s","setup","Choose setup manually by name",false,"", &allowedsetupnames);
    auto cmd_setupOptions = cmd.add<TCLAP::MultiArg<string>>("S","setup_options","Options for setup, key=value",false,"");

    auto cmd_maxevents = cmd.add<TCLAP::MultiArg<int>>("m","maxevents","Process only max events",false,"maxevents");

    TCLAP::ValuesConstraintExtra<decltype(analysis::PhysicsRegistry::GetList())> allowedPhysics(analysis::PhysicsRegistry::GetList());
    auto cmd_physicsclasses  = cmd.add<TCLAP::MultiArg<string>>("p","physics","Physics class to run", false, &allowedPhysics);

    auto cmd_output = cmd.add<TCLAP::ValueArg<string>>("o","output","Output file",false,"","filename");

    auto cmd_physicsOptions = cmd.add<TCLAP::MultiArg<string>>("O","options","Options for all physics classes, key=value",false,"");
    auto cmd_physicsclasses_opt = cmd.add<TCLAP::MultiArg<string>>("P","physics-opt","Physics class to run, with options: PhysicsClass:key=val,key=val", false, "");

    auto cmd_batchmode = cmd.add<TCLAP::MultiSwitchArg>("b","batch","Run in batch mode (no ROOT shell afterwards)",false);

    auto cmd_calibrations  = cmd.add<TCLAP::MultiArg<string>>("c","calibration","Calibration to run",false,"calibration");

    auto cmd_u_disablerecon  = cmd.add<TCLAP::SwitchArg>("","u_disablereconstruct","Unpacker: Disable Reconstruct (disables also all analysis)",false);

    auto cmd_p_disableParticleID  = cmd.add<TCLAP::SwitchArg>("","p_disableParticleID","Physics: Disable ParticleID",false);
    auto cmd_p_simpleParticleID  = cmd.add<TCLAP::SwitchArg>("","p_simpleParticleID","Physics: Use simple ParticleID (just protons/photons)",false);



    cmd.parse(argc, argv);
    if(cmd_verbose->isSet()) {
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());
    }

    // progress updates only when running interactively
    if(std_ext::system::isInteractive())
        ProgressCounter::Interval = 3;

    // enable caching of the calibration database
    ant::calibration::DataBase::OnDiskLayout::EnableCaching = true;

    // check if input files are readable
    for(const auto& inputfile : cmd_input->getValue()) {
        string errmsg;
        if(!std_ext::system::testopen(inputfile, errmsg)) {
            LOG(ERROR) << "Cannot open inputfile '" << inputfile << "': " << errmsg;
            return EXIT_FAILURE;
        }
    }

    // parse the setup options and tell the registry
    std::shared_ptr<OptionsList> setup_opts = make_shared<OptionsList>();
    if(cmd_setupOptions->isSet()) {
        for(const auto& opt : cmd_setupOptions->getValue()) {
            setup_opts->SetOption(opt);
        }
        ant::expconfig::SetupRegistry::SetSetupOptions(setup_opts);
    }


    // build the list of ROOT files first
    auto rootfiles = make_shared<WrapTFileInput>();
    for(const auto& inputfile : cmd_input->getValue()) {
        VLOG(5) << "ROOT File Manager: Looking at file " << inputfile;
        try {
            rootfiles->OpenFile(inputfile);
        } catch (const WrapTFile::ENotARootFile& e) {
            VLOG(5) <<  e.what();
        }
    }

    // check if there's a previous AntHeader present,
    // which could tell us the SetupName
    if(!cmd_setup->isSet()) {
        TAntHeader* previous_AntHeader;
        if(rootfiles->GetObject<TAntHeader>("AntHeader",previous_AntHeader)) {
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


    // now we can try to open the files with an unpacker
    std::unique_ptr<Unpacker::Module> unpacker = nullptr;
    for(const auto& inputfile : cmd_input->getValue()) {
        VLOG(5) << "Unpacker: Looking at file " << inputfile;
        try {
            auto unpacker_ = Unpacker::Get(inputfile);
            if(unpacker != nullptr && unpacker_ != nullptr) {
                LOG(ERROR) << "Can only handle one unpacker, but given input files suggest to use more than one.";
                return EXIT_FAILURE;
            }
            LOG(INFO) << "Found unpacker for file " << inputfile;
            unpacker = move(unpacker_);
        }
        catch(Unpacker::Exception& e) {
            VLOG(5) << "Unpacker: " << e.what();
        }
        catch(RawFileReader::Exception& e) {
            LOG(WARNING) << "RawFileReader: Error opening file " << inputfile << ": " << e.what();
        }
        catch(ExpConfig::ExceptionNoSetup&) {
            LOG(ERROR) << "The inputfile " << inputfile << " cannot be unpacked without a manually specified setupname. "
                       << "Consider using " << cmd_setup->longID();
            return EXIT_FAILURE;
        }
    }


    // we can finally we can create the available input readers
    // for the analysis

    list< unique_ptr<analysis::input::DataReader> > readers;

    // turn the unpacker into a input::DataReader
    {
        std::unique_ptr<Reconstruct_traits> reconstruct;
        if(!cmd_u_disablerecon->isSet()) {
            try {
                reconstruct = std_ext::make_unique<Reconstruct>();
            }
            catch(ExpConfig::ExceptionNoSetup&) {
                LOG(WARNING) << "Cannot activate reconstruct without setup";
            }
        }
        readers.push_back(std_ext::make_unique<analysis::input::AntReader>(
                              rootfiles,
                              move(unpacker),
                              move(reconstruct)
                              )
                          );
    }
    readers.push_back(std_ext::make_unique<analysis::input::PlutoReader>(rootfiles));
    readers.push_back(std_ext::make_unique<analysis::input::GoatReader>(rootfiles));


    // create the list of enabled calibrations here,
    // because now the readers (and underlying unpackers) did the work
    // of finding the config, so that we can simply ask the ExpConfig now
    list<shared_ptr<Calibration::PhysicsModule>> enabled_calibrations;
    if(cmd_calibrations->isSet()) {
        auto& setup = ExpConfig::Setup::Get();
        stringstream ss_calibrations;

        std::vector<std::string> calibrationnames = cmd_calibrations->getValue();
        std::list<std::string> leftovers(calibrationnames.begin(), calibrationnames.end());
        for(const auto& calibration : setup.GetCalibrations()) {
            ss_calibrations << calibration->GetName() << " ";
            if(!std_ext::contains(calibrationnames, calibration->GetName())) {
                continue;
            }
            leftovers.remove(calibration->GetName());
            enabled_calibrations.emplace_back(move(calibration));
        }

        if(!leftovers.empty()) {
            LOG(INFO) << "Available calibrations: " << ss_calibrations.str();
            for(string leftover : leftovers)
                LOG(ERROR) << "Specified calibration '" << leftover << "' not found in list of available calibrations";
            return EXIT_FAILURE;
        }
    }

    if( enabled_calibrations.size()
            +cmd_physicsclasses->getValue().size()
            +cmd_physicsclasses_opt->getValue().size() == 0) {
        LOG(ERROR) << "Please specify at least one physics module or enable at least one calibration";
        return EXIT_FAILURE;
    }

    // the real output file, create it here to get all
    // further ROOT objects into this output file
    unique_ptr<WrapTFileOutput> masterFile;
    if(cmd_output->isSet()) {
        // cd into masterFile upon creation
        masterFile = std_ext::make_unique<WrapTFileOutput>(cmd_output->getValue(), true);
    }

    // add the physics/calibrationphysics modules
    analysis::PhysicsManager pm(addressof(interrupt));
    std::shared_ptr<OptionsList> popts = make_shared<OptionsList>();

    if(cmd_physicsOptions->isSet()) {
        for(const auto& opt : cmd_physicsOptions->getValue()) {
            popts->SetOption(opt);
        }
    }

    for(const auto& classname : cmd_physicsclasses->getValue()) {
        try {
            auto physicsclass = analysis::PhysicsRegistry::Create(classname, popts);
            const auto instancename = physicsclass->GetName();
            pm.AddPhysics( move(physicsclass) );
            LOG(INFO) << "Activated physics class '" << classname << "' as instance '" << instancename << "'";
        } catch (const std::exception& e) {
            LOG(ERROR) << "Error while activating physics class \"" << classname << "\": " << e.what();
            return EXIT_FAILURE;
        } catch (...) {
            LOG(ERROR) << "Could not activate physics class for unknown reason";
            return EXIT_FAILURE;
        }
    }

    for(const auto& line : cmd_physicsclasses_opt->getValue()) {

        auto colpos = line.find(":");

        auto physicsname = line.substr(0,colpos);
        auto optstr = line.substr(colpos+1,line.npos);
        auto options = make_shared<OptionsList>(popts);
        options->SetOptions(optstr);
        try {
            auto physicsclass = analysis::PhysicsRegistry::Create(physicsname, options);
            const auto instancename = physicsclass->GetName();
            pm.AddPhysics( move(physicsclass) );
            LOG(INFO) << "Activated physics class '" << physicsname << "' as instance '" << instancename << "' with options " << optstr;
        } catch (...) {
            LOG(ERROR) << "Physics class '" << line << "' not found";
            return EXIT_FAILURE;
        }

        auto unused_popts = options->GetUnused();
        if(!unused_popts.empty()) {
            LOG(ERROR) << "Physics class '" << physicsname << "' did not use the options: " << unused_popts;
            LOG(INFO)  << "Did you mean: " << options->GetNotFound();
            return EXIT_FAILURE;
        }
    }

    for(const auto& calibration : enabled_calibrations) {
        const auto& physicsclasses = calibration->GetPhysicsModules();
        if(physicsclasses.empty()) {
            // this is actually more an implementation error...
            LOG(ERROR) << "Calibration '" << calibration->GetName() << "' did not specify any physics classes";
            return EXIT_FAILURE;
        }
        for(const std::string classname : physicsclasses) {
            try {
                pm.AddPhysics( analysis::PhysicsRegistry::Create(classname, popts) );
                LOG(INFO) << "Activated physics class '" << classname << "'";
            } catch (...) {
                LOG(ERROR) << "Physics class '" << classname << "' requested by calibration '"
                           << calibration->GetName()
                           << "' not found";
                return EXIT_FAILURE;
            }
        }
    }

    // check global unused after activating all physics classes
    auto global_unused_popts =  popts->GetUnused();
    if(!global_unused_popts.empty()) {
        LOG(ERROR) << "The following global physics options were not used by any activated physics class: " << global_unused_popts;
        LOG(INFO)  << "Did you mean: " << popts->GetNotFound();
        return EXIT_FAILURE;
    }

    // check setup options after physics initialization
    auto unused_setup_opts = setup_opts->GetUnused();
    if(!unused_setup_opts.empty()) {
        LOG(ERROR) << "The following setup options were not used: " << unused_setup_opts;
        LOG(INFO)  << "Did you mean: " << setup_opts->GetNotFound();
        return EXIT_FAILURE;
    }

    // set up particle ID

    if(!cmd_p_disableParticleID->isSet()) {
        try {
            unique_ptr<analysis::utils::ParticleID> particleID;
            if(cmd_p_simpleParticleID->isSet()) {
                particleID = std_ext::make_unique<analysis::utils::SimpleParticleID>();
            } else {
                auto& setup = ExpConfig::Setup::Get();
                particleID = std_ext::make_unique<analysis::utils::CBTAPSBasicParticleID>(setup.GetPIDCutsDirectory());
            }
            analysis::utils::ParticleID::SetDefault(move(particleID));
        }
        catch(ExpConfig::ExceptionNoSetup&) {
            LOG(WARNING) << "Cannot activate particle ID without setup";
        }
    } else {
        LOG(INFO) << "ParticleID disabled by command line";
    }

    // create some variables for running
    long long maxevents = cmd_maxevents->isSet()
            ? cmd_maxevents->getValue().back()
            :  numeric_limits<long long>::max();


    // this method does the hard work...
    pm.ReadFrom(move(readers), maxevents);
    rootfiles = nullptr; // cleanup opened ROOT files for reading

    TAntHeader* header = new TAntHeader();
    gDirectory->Add(header);
    {
        const auto& tidRange = pm.GetProcessedTIDRange();
        header->FirstID = tidRange.Start();
        header->LastID  = tidRange.Stop();
    }

    // add some more info about the current state
    try
    {
        auto& setup = ExpConfig::Setup::Get();
        header->SetupName = setup.GetName();
        if(auto calmgr = setup.GetCalibrationDataManager()) {
            GitInfo gitinfo_db(calmgr->GetCalibrationDataFolder());
            header->GitInfoDatabase = gitinfo_db.GetDescription();
        }
    }
    catch(ExpConfig::ExceptionNoSetup&) {
        LOG(WARNING) << "HeaderInfo not complete without setup";
    }

    GitInfo gitinfo;
    header->GitInfo = gitinfo.GetDescription();
    header->WorkingDir = std_ext::system::getCwd();
    header->CmdLine = std_ext::system::buildCmdLine(argc, argv);
    VLOG(5) << "Added header: " << *header;

    if(terminated)
        return EXIT_FAILURE+1;

    if(!cmd_batchmode->isSet()) {
        if(!std_ext::system::isInteractive()) {
            LOG(INFO) << "No TTY attached. Not starting ROOT shell.";
        }
        else {
            argc=0; // prevent TRint to parse any cmdline
            auto app = new TRint("Ant",&argc,argv,nullptr,0,true);

            if(masterFile != nullptr)
                LOG(INFO) << "Stopped running, but close ROOT properly to write data to disk.";

            pm.ShowResults();
            app->Run(kTRUE); // really important to return...
            masterFile = nullptr;   // and to destroy the master WrapTFile before TRint is destroyed
            // call this before application tear down
            gROOT->EndOfProcessCleanups();
            // delete app; // do not delete, fixes segfault on some machines
        }

    }

    return EXIT_SUCCESS;
}
