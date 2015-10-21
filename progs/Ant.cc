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

#include "unpacker/Unpacker.h"
#include "unpacker/RawFileReader.h"

#include "reconstruct/Reconstruct.h"

#include "tree/UnpackerReader.h"
#include "tree/UnpackerWriter.h"
#include "tree/THeaderInfo.h"
#include "tree/TAntHeader.h"

#include "base/std_ext/vector.h"
#include "base/WrapTFile.h"
#include "base/Logger.h"
#include "base/CmdLine.h"
#include "base/detail/tclap/ValuesConstraintExtra.h"
#include "base/WrapTFile.h"
#include "base/std_ext/system.h"
#include "base/GitInfo.h"


#include "TRint.h"
#include "TSystem.h"

#include <sstream>
#include <string>

using namespace std;
using namespace ant;

bool running = true;
class MyTInterruptHandler : public TSignalHandler {
public:
    MyTInterruptHandler() : TSignalHandler(kSigInterrupt, kFALSE) { }

    Bool_t  Notify() {
        if (fDelay) {
            fDelay++;
            return kTRUE;
        }
        running = false;
        cout << " >>> Interrupted! " << endl;
        return kTRUE;
    }
};

int main(int argc, char** argv) {
    SetupLogger();

    // check for bash completion commands
    if(argc >= 2) {
        const string arg1(argv[1]);

        if(arg1 == "--list-physics") {
            for(const auto& name : analysis::PhysicsRegistry::GetList()) {
                cout << name << endl;
            }
            return 0;
        }

        if(arg1 == "--list-setups") {
            for(const auto& name : ExpConfig::Setup::GetNames()) {
                cout << name << endl;
            }
            return 0;
        }

        if(arg1 == "--list-calibrations") {

            if(argc == 3) {

                const string setup_name(argv[2]);

                ExpConfig::Setup::ManualName = setup_name;
                const auto setup = ExpConfig::Setup::GetLastFound();

                if(setup) {
                    for(const auto& calibration : setup->GetCalibrations()) {
                        cout << calibration->GetName() << endl;
                    }
                }
            }

            return 0;
        }
    }

    TCLAP::CmdLine cmd("Ant", ' ', "0.1");

    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"int");
    auto cmd_input  = cmd.add<TCLAP::MultiArg<string>>("i","input","Input files",true,"filename");

    TCLAP::ValuesConstraintExtra<decltype(ExpConfig::Setup::GetNames())> allowedsetupnames(ExpConfig::Setup::GetNames());
    auto cmd_setup  = cmd.add<TCLAP::ValueArg<string>>("s","setup","Choose setup manually by name",false,"", &allowedsetupnames);
    auto cmd_setupOptions = cmd.add<TCLAP::MultiArg<string>>("S","setup_options","Options for setup, key=value",false,"");

    auto cmd_maxevents = cmd.add<TCLAP::ValueArg<int>>("m","maxevents","Process only max events",false, 0, "maxevents");

    TCLAP::ValuesConstraintExtra<decltype(analysis::PhysicsRegistry::GetList())> allowedPhysics(analysis::PhysicsRegistry::GetList());
    auto cmd_physicsclasses  = cmd.add<TCLAP::MultiArg<string>>("p","physics","Physics class to run", false, &allowedPhysics);

    auto cmd_output = cmd.add<TCLAP::ValueArg<string>>("o","output","Output file",false,"","filename");

    auto cmd_physicsOptions = cmd.add<TCLAP::MultiArg<string>>("O","options","Options for all physics classes, key=value",false,"");
    auto cmd_physicsclasses_opt = cmd.add<TCLAP::MultiArg<string>>("P","physics-opt","Physics class to run, with options: PhysicsClass:key=val,key=val", false, "");

    auto cmd_batchmode = cmd.add<TCLAP::SwitchArg>("b","batch","Run in batch mode (no ROOT shell afterwards)",false);

    auto cmd_calibrations  = cmd.add<TCLAP::MultiArg<string>>("c","calibration","Calibration to run",false,"calibration");

    auto cmd_unpackerout  = cmd.add<TCLAP::ValueArg<string>>("u","unpackerout","Unpacker output file",false,"","outputfile");
    auto cmd_u_writeuncal  = cmd.add<TCLAP::SwitchArg>("","u_writeuncalibrated","Unpacker: Output UNcalibrated detector reads (before reconstruct)",false);
    auto cmd_u_disablerecon  = cmd.add<TCLAP::SwitchArg>("","u_disablereconstruct","Unpacker: Disable Reconstruct (disables also all analysis)",false);
    auto cmd_u_writecal  = cmd.add<TCLAP::SwitchArg>("","u_writecalibrated","Unpacker: Output calibrated detector reads (only if Reconstruct found)",false);

    auto cmd_p_disableParticleID  = cmd.add<TCLAP::SwitchArg>("","p_disableParticleID","Physics: Disable ParticleID",false);
    auto cmd_p_simpleParticleID  = cmd.add<TCLAP::SwitchArg>("","p_simpleParticleID","Physics: Use simple ParticleID (just protons/photons)",false);



    cmd.parse(argc, argv);
    if(cmd_verbose->isSet()) {
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());
    }

    if(std_ext::system::isInteractive())
        RawFileReader::OutputPerformanceStats = 3;




    int fake_argc=0;
    char** fake_argv=nullptr;
    TRint app("Ant",&fake_argc,fake_argv,nullptr,0,true);
    auto oldsig = app.GetSignalHandler();
    oldsig->Remove();
    auto mysig = new MyTInterruptHandler();
    mysig->Add();
    gSystem->AddSignalHandler(mysig);

    // check if input files are readable
    for(const auto& inputfile : cmd_input->getValue()) {
        string errmsg;
        if(!std_ext::system::testopen(inputfile, errmsg)) {
            LOG(ERROR) << "Cannot open inputfile '" << inputfile << "': " << errmsg;
            return 1;
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
        } catch (const std::runtime_error& e) {
            VLOG(5) << "Could not open ROOT file " << inputfile;
        }
    }

    // then init the unpacker root input file manager
    auto unpackerFile = std_ext::make_unique<tree::UnpackerReader>(rootfiles);

    // search for header info?
    if(unpackerFile->OpenInput()) {
        LOG(INFO) << "Found complete set of input ROOT trees for unpacker";
        THeaderInfo headerInfo;
        if(unpackerFile->GetUniqueHeaderInfo(headerInfo)) {
            VLOG(5) << "Found unique header info " << headerInfo;
            if(!headerInfo.SetupName.empty()) {
                ExpConfig::Setup::ManualName = headerInfo.SetupName;
                LOG(INFO) << "Using header info to manually set the setup name to " << ExpConfig::Setup::ManualName;
            }
        }
    }

    // override the setup name from cmd line
    if(cmd_setup->isSet()) {
        ExpConfig::Setup::ManualName = cmd_setup->getValue();
        if(ExpConfig::Setup::ManualName.empty())
            LOG(INFO) << "Commandline override to auto-search for setup config (might fail)";
        else
            LOG(INFO) << "Commandline override setup name to '" << ExpConfig::Setup::ManualName << "'";
    }


    // now we can try to open the files with an unpacker
    std::unique_ptr<Unpacker::Reader> unpacker = nullptr;
    for(const auto& inputfile : cmd_input->getValue()) {
        VLOG(5) << "Unpacker: Looking at file " << inputfile;
        try {
            auto unpacker_ = Unpacker::Get(inputfile);
            if(unpacker != nullptr && unpacker_ != nullptr) {
                LOG(ERROR) << "Can only handle one unpacker, but given input files suggest to use more than one.";
                return 1;
            }
            LOG(INFO) << "Found unpacker for file " << inputfile;
            unpacker = move(unpacker_);
        }
        catch(Unpacker::Exception e) {
            VLOG(5) << "Unpacker: " << e.what();
        }
        catch(RawFileReader::Exception e) {
            LOG(WARNING) << "Unpacker: Error opening file "<<inputfile<<": " << e.what();
        }
        catch(ExpConfig::ExceptionNoConfig) {
            LOG(ERROR) << "The inputfile " << inputfile << " cannot be unpacked without a manually specified setupname. "
                       << "Consider using " << cmd_setup->longID();
            return 1;
        }
    }

    // select the right source for the unpacker stage
    if(unpacker != nullptr && unpackerFile->IsOpen()) {
        LOG(WARNING) << "Found file suitable for unpacker and ROOT file for unpacker stage, preferring raw data file";
    }
    else if(unpackerFile->IsOpen()) {
        LOG(INFO) << "Running unpacker stage from input ROOT file(s)";
        unpacker = move(unpackerFile);
    }
    else if(unpacker != nullptr) {
        LOG(INFO) << "Running unpacker stage from raw data file";
    }


    // we can finally we can create the available input readers
    // for the analysis

    list< unique_ptr<analysis::input::DataReader> > readers;

    if(unpacker) {
        // turn the unpacker into a input::DataReader
        auto reconstruct = cmd_u_disablerecon->isSet() ? nullptr : std_ext::make_unique<Reconstruct>();
        auto unpacker_reader =
                std_ext::make_unique<analysis::input::AntReader>(
                    move(unpacker),
                    move(reconstruct)
                    );
        // it may write stuff during the unpacker stage
        if(cmd_unpackerout->isSet()) {
            unpacker_reader->EnableUnpackerWriter(
                        cmd_unpackerout->getValue(),
                        cmd_u_writeuncal->isSet(),
                        cmd_u_writecal->isSet()
                        );
        }
        readers.push_back(move(unpacker_reader));
    }

    // the unpackers might have figured out what setup to use...
    auto setup = ExpConfig::Setup::GetLastFound();
    auto tagger = setup ? setup->GetDetector<TaggerDetector_t>() : nullptr;

    readers.push_back(std_ext::make_unique<analysis::input::PlutoReader>(rootfiles, tagger));
    readers.push_back(std_ext::make_unique<analysis::input::GoatReader>(rootfiles));

    // create the list of enabled calibrations here,
    // because now the readers (and underlying unpackers) did the work
    // of finding the config, so that we can simply ask the ExpConfig now
    list<shared_ptr<Calibration::PhysicsModule>> enabled_calibrations;
    if(cmd_calibrations->isSet()) {
        if(setup==nullptr) {
            stringstream ss_setups;
            for(auto name : ExpConfig::Setup::GetNames()) {
                ss_setups << name << " ";
            }
            LOG(INFO)  << "Available setups: " << ss_setups.str();
            LOG(ERROR) << "Please specify a --setup if you want to use calibrations as physics modules";
            return 1;
        }
        else {
            stringstream ss_calibrations;

            std::vector<std::string> calibrationnames = cmd_calibrations->getValue();
            std::list<std::string> leftovers(calibrationnames.begin(), calibrationnames.end());
            for(const auto& calibration : setup->GetCalibrations()) {
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
                return 1;
            }
        }
    }

    if(enabled_calibrations.size()+cmd_physicsclasses->getValue().size() == 0) {
        LOG(ERROR) << "Please specify at least one physics module or enable at least one calibration";
        return 1;
    }

    // the real output file, create it here to get all
    // further ROOT objects into this output file
    unique_ptr<WrapTFileOutput> masterFile;
    if(cmd_output->isSet()) {
        masterFile = std_ext::make_unique<WrapTFileOutput>(cmd_output->getValue(),
                                                    WrapTFileOutput::mode_t::recreate,
                                                     true); // cd into masterFile upon creation
    }

    // add the physics/calibrationphysics modules
    analysis::PhysicsManager pm;
    std::shared_ptr<OptionsList> popts = make_shared<OptionsList>();

    if(cmd_physicsOptions->isSet()) {
        for(const auto& opt : cmd_physicsOptions->getValue()) {
            popts->SetOption(opt);
        }
    }

    for(const auto& classname : cmd_physicsclasses->getValue()) {
        try {
            pm.AddPhysics( analysis::PhysicsRegistry::Create(classname, popts) );
            LOG(INFO) << "Activated physics class '" << classname << "'";
        } catch (...) {
            LOG(ERROR) << "Physics class '" << classname << "' is not found.";
            return 1;
        }
    }

    for(const auto& line : cmd_physicsclasses_opt->getValue()) {

        auto colpos = line.find(":");

        auto classname = line.substr(0,colpos);
        auto optstr = line.substr(colpos+1,line.npos);
        auto options = make_shared<OptionsList>(popts);
        options->SetOptions(optstr);
        try {
            pm.AddPhysics( analysis::PhysicsRegistry::Create(classname, options) );
            LOG(INFO) << "Activated physics class '" << classname << "'";
        } catch (...) {
            LOG(ERROR) << "Physics class '" << line << "' is not found.";
            return 1;
        }
    }


    for(const auto& calibration : enabled_calibrations) {
        pm.AddPhysics(calibration->GetPhysicsModule());
    }



    // set up particle ID

    if(!cmd_p_disableParticleID->isSet()) {

        unique_ptr<analysis::utils::ParticleID> particleID;

        if(cmd_p_simpleParticleID->isSet())
        {
            particleID = std_ext::make_unique<analysis::utils::SimpleParticleID>();
        }
        else
        {
            if(auto setup = ExpConfig::Setup::GetLastFound()) {
                particleID = std_ext::make_unique<analysis::utils::CBTAPSBasicParticleID>
                             (setup->GetPIDCutsDirectory());
            } else {
                LOG(WARNING) << "No Setup found while loading ParticleID cuts.";
            }
        }
        pm.SetParticleID(move(particleID));

    } else {
        LOG(INFO) << "ParticleID disabled by command line";
    }



    // create some variables for running
    long long maxevents = cmd_maxevents->isSet()
            ? cmd_maxevents->getValue()
            :  numeric_limits<long long>::max();
    TAntHeader* header = new TAntHeader();

    // this method does the hard work...
    pm.ReadFrom(move(readers), maxevents, running, *header);

    // add some more info about the current state
    if(auto setup = ExpConfig::Setup::GetLastFound()) {
        header->SetupName = setup->GetName();
    }
    header->GitInfo = GitInfo::GetDescription();
    if(!header->GitInfo.empty()) {
        VLOG(5) << "Added git info: " << header->GitInfo;
    }
    char* pwd = get_current_dir_name();
    header->WorkingDir = pwd;
    free(pwd);
    VLOG(5) << "Added working directory: " << header->WorkingDir;
    while(--argc>0) {
        header->CmdLine += *argv++;
        header->CmdLine += " ";
    }
    VLOG(5) << "Added command line: " << header->CmdLine;

    if(!cmd_batchmode->isSet()) {
        if(!std_ext::system::isInteractive()) {
            LOG(INFO) << "No TTY attached. Not starting ROOT shell.";
        } else {

            mysig->Remove();
            oldsig->Add();
            gSystem->AddSignalHandler(oldsig);
            delete mysig;

            if(masterFile != nullptr)
                LOG(INFO) << "Stopped running, but close ROOT properly to write data to disk.";

            pm.ShowResults();
            app.Run(kTRUE); // really important to return...
            masterFile = nullptr;   // and to destroy the master WrapTFile before TRint is destroyed
        }

    }

    return 0;
}
