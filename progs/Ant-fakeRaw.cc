#include "expconfig/ExpConfig.h"
#include "analysis/input/ant/AntReader.h"
#include "unpacker/UnpackerAcqu.h"

#include "expconfig/setups/Setup.h"

#include "tree/TEvent.h"
#include "tree/TEventData.h"
#include "tree/TAntHeader.h"

#include "base/WrapTFile.h"
#include "base/std_ext/system.h"
#include "base/CmdLine.h"
#include "base/Logger.h"

#include <memory>
#include <signal.h>

using namespace std;
using namespace ant;

static volatile bool interrupt = false;

int main(int argc, char** argv) {
    SetupLogger();

    signal(SIGINT, [] (int) {
        interrupt = true;
    });

    TCLAP::CmdLine cmd("Ant-fakeRaw", ' ', "0.1");

    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"int");
    auto cmd_setup  = cmd.add<TCLAP::ValueArg<string>>("s","setup","Choose setup manually by name",false,"","setup");
    auto cmd_input  = cmd.add<TCLAP::ValueArg<string>>("i","input","Input files",true,"","filename");
    auto cmd_output = cmd.add<TCLAP::ValueArg<string>>("o","output","Output file",false,"","filename");

    cmd.parse(argc, argv);
    if(cmd_verbose->isSet()) {
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());
    }

    // check if input file is readable
    const string& inputfile = cmd_input->getValue();
    {
        string errmsg;
        if(!std_ext::system::testopen(inputfile, errmsg)) {
            LOG(ERROR) << "Cannot open inputfile '" << inputfile << "': " << errmsg;
            return 1;
        }
    }

    auto inputrootfile = make_shared<WrapTFileInput>(inputfile);

    // check if there's a previous AntHeader present,
    // which could tell us the SetupName
    TAntHeader* previous_AntHeader = nullptr;
    if(inputrootfile->GetObject<TAntHeader>("AntHeader",previous_AntHeader)) {
        const auto& setupname = previous_AntHeader->SetupName;
        if(!setupname.empty()) {
            ExpConfig::Setup::SetManualName(setupname);
            LOG(INFO) << "Setup name set to '" << setupname << "' from input file";
        }
        else
            LOG(WARNING) << "Found AntHeader in input files, but SetupName was empty";
    }

    // override the setup name from cmd line
    if(cmd_setup->isSet()) {
        const auto& setupname = cmd_setup->getValue();
        ExpConfig::Setup::SetManualName(setupname, false);
        LOG(INFO) << "Commandline override setup name to '" << setupname << "'";
    }


    analysis::input::AntReader reader(inputrootfile, nullptr, nullptr);


    unsigned nEvents = 0;
    TEvent event;
    while(reader.ReadNextEvent(event)) {
        if(interrupt)
            break;
        LOG(INFO) << event.Reconstructed().DetectorReadHits.size();
        nEvents++;
    }
    LOG(INFO) << nEvents << " events processed";



    return EXIT_SUCCESS;
}
