
#include "analysis/input/goat/GoatReader.h"
#include "analysis/OutputManager.h"
#include "analysis/physics/Physics.h"
#include "analysis/physics/omega/omega.h"
#include "analysis/physics/common/DataOverview.h"

#include "base/Logger.h"
#include "base/CmdLine.h"
#include "base/ReadTFiles.h"
#include "TRint.h"

using namespace std;
using namespace ant::output;
using namespace ant;
using namespace ant::analysis;

int main(int argc, char** argv) {
    SetupLogger();


    TCLAP::CmdLine cmd("Omega Analysis", ' ', "0.1");
    auto input  = cmd.add<TCLAP::MultiArg<string>>("i","input","GoAT input files",true,"string");
    auto output = cmd.add<TCLAP::ValueArg<string>>("o","output","Output file",false,"","string");
    auto max_event = cmd.add<TCLAP::ValueArg<int>>("","stop-at","Stop at event number",false,0,"int");
    auto batchmode = cmd.add<TCLAP::SwitchArg>("b","batch","Run in batch mode (No ROOT Windows)",false);
    auto verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"int");

    cmd.parse(argc, argv);
    if(verbose->isSet()) {
        el::Loggers::setVerboseLevel(verbose->getValue());
    }


    OutputManager om;

    if(output->isSet())
        om.SetNewOutput(output->getValue());

    PhysicsManager pm;

    pm.AddPhysics<OmegaEtaG>(OmegaBase::DataMode::Reconstructed);
    pm.AddPhysics<DataOverview>();

    auto filemanager = make_shared<ReadTFiles>();

    for(const auto& file : input->getValue())
        filemanager->OpenFile(file);


    input::GoatReader reader(filemanager);



//    long long maxevents = max_event->isSet()
//            ? max_event->getValue()
//            :  numeric_limits<long long>::max();
//    // this method does the hard work...
//    bool running = true;
//    pm.ReadFrom(move(readers), maxevents, running);

    if(!batchmode->isSet()) {
        int a=0;
        char** b=nullptr;
        TRint app("omega",&a,b);
        pm.ShowResults();
        app.Run(kTRUE);
    }

    return 0;

}

