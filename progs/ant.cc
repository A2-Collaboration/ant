

#include "analysis/input/DataReader.h"
#include "analysis/input/ant/AntReader.h"
#include "analysis/input/goat/GoatReader.h"
#include "analysis/OutputManager.h"

#include "analysis/physics/Physics.h"
#include "analysis/physics/omega/omega.h"
#include "analysis/physics/common/DataOverview.h"
#include "analysis/physics/common/CandidatesAnalysis.h"

#include "expconfig/ExpConfig.h"

#include "unpacker/Unpacker.h"
#include "unpacker/RawFileReader.h"

#include "tree/UnpackerReader.h"
#include "tree/THeaderInfo.h"

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

using namespace std;
using namespace ant::output;
using namespace ant;
using namespace ant::analysis;





int main(int argc, char** argv) {
    SetupLogger();

    TCLAP::CmdLine cmd("ant", ' ', "0.1");
    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"int");
    auto cmd_input  = cmd.add<TCLAP::MultiArg<string>>("i","input","Input files",true,"string");
    auto cmd_setup  = cmd.add<TCLAP::ValueArg<string>>("s","setup","Choose setup",false,"","string");



    cmd.parse(argc, argv);
    if(cmd_verbose->isSet()) {
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());
    }

    // build the general ROOT file manager first
    auto filemanager = make_shared<ReadTFiles>();
    for(const auto& inputfile : cmd_input->getValue()) {
        VLOG(5) << "ROOT File Manager: Looking at file " << inputfile;
        if(filemanager->OpenFile(inputfile))
            LOG(INFO) << "Opened file '" << inputfile << "' as ROOT file";
        else
            VLOG(5) << "Could not add " << inputfile << " to ROOT file manager";
    }

    // then init the unpacker root input file manager
    auto unpackerFile = std_ext::make_unique<tree::UnpackerReader>(filemanager);

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
            VLOG(5) << "Unpacker::Get said: " << e.what();
            unpacker = nullptr;
        }
        catch(RawFileReader::Exception e) {
            LOG(WARNING) << "Error opening file "<<inputfile<<": " << e.what();
            unpacker = nullptr;
        }
        catch(ExpConfig::ExceptionNoConfig) {
            LOG(ERROR) << "The file "<<inputfile<<" cannot be unpacked without a manually specified setupname. "
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


    while(auto item = unpacker->NextItem()) {
        cout << *item << endl;
    }



    //unpackerFile.DetectorRead->Class()->

    //cout << unpackerFile.HeaderInfo->Class_Name() << endl;

    return 0;

}

