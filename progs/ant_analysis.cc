
#include "analysis/input/ant/AntReader.h"
#include "analysis/OutputManager.h"

#include "analysis/physics/Physics.h"
#include "analysis/physics/omega/omega.h"
#include "analysis/physics/common/DataOverview.h"
#include "analysis/physics/common/TrackAnalysis.h"

#include "base/Logger.h"
#include "base/detail/CmdLine.h"

#include "TRint.h"

#include <string>

using namespace std;
using namespace ant::output;
using namespace ant;
using namespace ant::analysis;

int main(int argc, char** argv) {
    SetupLogger();


    TCLAP::CmdLine cmd("Omega Analysis", ' ', "0.1");
    auto input  = cmd.add<TCLAP::MultiArg<string>>("i","input","ant root file (events)",true,"string");
    auto physicsclasses  = cmd.add<TCLAP::MultiArg<string>>("p","physics","Physics Class to run",true,"string");
    auto output = cmd.add<TCLAP::ValueArg<string>>("o","output","Output file",false,"","string");
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

    for( auto& p : physicsclasses->getValue()) {
        if(p == "OmetaEtaG")
            pm.AddPhysics<OmegaEtaG>(OmegaBase::DataMode::Reconstructed);
        else if (p== "DataOverview")
            pm.AddPhysics<DataOverview>();
        else if (p=="TrackAnalysis")
            pm.AddPhysics<TrackAnalysis>();
        else {
            LOG(WARNING) << "Unknown physics class " << p;
        }
    }

    input::AntReader reader;

    for(auto& file : input->getValue())
            reader.AddInputFile(file);

    reader.Initialize();


    pm.ReadFrom(reader);

    if(!batchmode->isSet()) {
        int a=0;
        char** b=nullptr;
        TRint app("omega",&a,b);
        pm.ShowResults();
        app.Run(kTRUE);
    }

    return 0;

}

