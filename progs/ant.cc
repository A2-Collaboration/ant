

#include "analysis/input/DataReader.h"
#include "analysis/input/ant/AntReader.h"
#include "analysis/input/goat/GoatReader.h"
#include "analysis/input/pluto/PlutoReader.h"
#include "analysis/OutputManager.h"

#include "analysis/physics/Physics.h"
#include "analysis/physics/omega/omega.h"
#include "analysis/physics/common/DataOverview.h"
#include "analysis/physics/common/CandidatesAnalysis.h"

#include "expconfig/ExpConfig.h"

#include "unpacker/Unpacker.h"
#include "unpacker/RawFileReader.h"

#include "reconstruct/Reconstruct.h"

#include "tree/UnpackerReader.h"
#include "tree/UnpackerWriter.h"
#include "tree/THeaderInfo.h"
#include "tree/TDetectorRead.h"
#include "tree/TEvent.h"

#include "base/std_ext.h"
#include "base/Logger.h"
#include "base/CmdLine.h"
#include "base/ReadTFiles.h"

#include "TRint.h"
#include "TClass.h"
#include "TTree.h"
#include "TFile.h"

#include <sstream>
#include <string>
#include <chrono>
#include <cstdio>
#include <cerrno>


using namespace std;
using namespace ant::output;
using namespace ant;
using namespace ant::analysis;


bool running = true;
void myCrashHandler(int sig);


int main(int argc, char** argv) {
    SetupLogger();
    el::Helpers::setCrashHandler(myCrashHandler);

    TCLAP::CmdLine cmd("ant", ' ', "0.1");
    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"level");
    auto cmd_input  = cmd.add<TCLAP::MultiArg<string>>("i","input","Input files",true,"inputfile");
    auto cmd_setup  = cmd.add<TCLAP::ValueArg<string>>("s","setup","Choose setup manually by name",false,"","setupname");
    auto cmd_maxevents = cmd.add<TCLAP::ValueArg<int>>("m","maxevents","Process only max events",false, 0, "maxevents");
    auto cmd_physicsclasses  = cmd.add<TCLAP::MultiArg<string>>("p","physics","Physics class to run",false,"physics");
    auto cmd_output = cmd.add<TCLAP::ValueArg<string>>("o","output","Output file",false,"","outputfile");
    auto cmd_batchmode = cmd.add<TCLAP::SwitchArg>("b","batch","Run in batch mode (no ROOT shell afterwards)",false);

    auto cmd_calibrations  = cmd.add<TCLAP::MultiArg<string>>("c","calibration","Calibration to run",false,"calibration");

    auto cmd_unpackerout  = cmd.add<TCLAP::ValueArg<string>>("u","unpackerout","Unpacker output file",false,"","outputfile");
    auto cmd_u_writeuncal  = cmd.add<TCLAP::SwitchArg>("","u_writeuncalibrated","Unpacker: Output UNcalibrated detector reads (before reconstruct)",false);
    auto cmd_u_disablerecon  = cmd.add<TCLAP::SwitchArg>("","u_disablereconstruct","Unpacker: Disable Reconstruct (disables also all analysis)",false);
    auto cmd_u_writecal  = cmd.add<TCLAP::SwitchArg>("","u_writecalibrated","Unpacker: Output calibrated detector reads (only if Reconstruct found)",false);



    cmd.parse(argc, argv);
    if(cmd_verbose->isSet()) {
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());
    }
    RawFileReader::OutputPerformanceStats = 3;

    // check if input files are readable
    for(const auto& inputfile : cmd_input->getValue()) {
        /// \todo find some non-C'ish way of doing this?
        FILE* fp = fopen(inputfile.c_str(),"r");
        if(fp==NULL) {
            LOG(ERROR) << "Inputfile '" << inputfile << "' could not be opened for reading: " << strerror(errno);
            return 1;
        }
        fclose(fp);
    }

    // build the list of ROOT files first
    auto rootfiles = make_shared<ReadTFiles>();
    for(const auto& inputfile : cmd_input->getValue()) {
        VLOG(5) << "ROOT File Manager: Looking at file " << inputfile;
        if(rootfiles->OpenFile(inputfile))
            LOG(INFO) << "Opened file '" << inputfile << "' as ROOT file";
        else
            VLOG(5) << "Could not add " << inputfile << " to ROOT file manager";
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
                ExpConfig::ManualSetupName = headerInfo.SetupName;
                LOG(INFO) << "Using header info to manually set the setup name to " << ExpConfig::ManualSetupName;
            }
        }
    }

    // override the setup name from cmd line
    if(cmd_setup->isSet()) {
        ExpConfig::ManualSetupName = cmd_setup->getValue();
        if(ExpConfig::ManualSetupName.empty())
            LOG(INFO) << "Commandline override to auto-search for setup config (might fail)";
        else
            LOG(INFO) << "Commandline override setup name to '" << ExpConfig::ManualSetupName << "'";
    }


    // now we can try to open the files with an unpacker
    std::unique_ptr<Unpacker::Reader> unpacker = nullptr;
    for(const auto& inputfile : cmd_input->getValue()) {
        VLOG(5) << "Unpacker: Looking at file " << inputfile;
        try {
            auto unpacker_ = Unpacker::Get(inputfile);
            if(unpacker != nullptr && unpacker_ != nullptr) {
                LOG(ERROR) << "Can only handle one unpacker, but input files suggest to use more than one.";
                return 1;
            }
            unpacker = move(unpacker_);
        }
        catch(Unpacker::Exception e) {
            VLOG(5) << "Unpacker: " << e.what();
            unpacker = nullptr;
        }
        catch(RawFileReader::Exception e) {
            LOG(WARNING) << "Unpacker: Error opening file "<<inputfile<<": " << e.what();
            unpacker = nullptr;
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

    list< unique_ptr<input::DataReader> > readers;

    if(unpacker) {
        // turn the unpacker into a input::DataReader
        auto reconstruct = cmd_u_disablerecon->isSet() ? nullptr : std_ext::make_unique<Reconstruct>();
        auto unpacker_reader =
                std_ext::make_unique<input::AntReader>(
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

    readers.push_back(std_ext::make_unique<input::PlutoReader>(rootfiles));
    readers.push_back(std_ext::make_unique<input::GoatReader>(rootfiles));

    // create the list of enabled calibrations here,
    // because now the readers (and underlying unpackers) did the work
    // of finding the config, so that
    list<shared_ptr<Calibration::PhysicsModule>> enabled_calibrations;
    if(cmd_calibrations->isSet()) {
        auto setup = ExpConfig::Setup::GetLastFound();
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

            for(string leftover : leftovers)
                LOG(WARNING) << "Specified calibration " << leftover << " not found in list of available calibrations";

            if(!leftovers.empty())
                LOG(INFO) << "Available calibrations: " << ss_calibrations.str();
        }
    }


    // the real output file, create it here to get all
    // further ROOT objects into this output file
    OutputManager om;
    if(cmd_output->isSet())
        om.SetNewOutput(cmd_output->getValue());

    PhysicsManager pm;
    for(const auto& classname : cmd_physicsclasses->getValue()) {
        try {
            pm.AddPhysics( PhysicsRegistry::Create(classname) );
            LOG(INFO) << "Activated physics class '" << classname << "'";
        } catch (...) {
            LOG(WARNING) << "Physics class '" << classname << "' is not found.";
        }
    }

    for(const auto& calibration : enabled_calibrations) {
        pm.AddPhysics(calibration->GetPhysicsModule());
    }

    long long maxevents = cmd_maxevents->isSet()
            ? cmd_maxevents->getValue()
            :  numeric_limits<long long>::max();
    // this method does the hard work...
    pm.ReadFrom(move(readers), maxevents, running);

    if(!cmd_batchmode->isSet()) {
        LOG(INFO) << "Stopped running, but close ROOT properly to write data to disk.";
        int a=0;
        char** b=nullptr;
        TRint app("ant",&a,b);
        pm.ShowResults();
        app.Run(kTRUE); // really important to return
    }

    return 0;
}

void myCrashHandler(int sig) {
    if(sig == SIGINT) {
        LOG(INFO) << "Ctrl-C received, stop running...";
        running = false;
        return;
    }
    // FOLLOWING LINE IS ABSOLUTELY NEEDED AT THE END IN ORDER TO ABORT APPLICATION
    el::Helpers::crashAbort(sig);
}