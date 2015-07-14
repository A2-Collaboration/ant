#include "input/goat/GoatReader.h"
#include "OutputManager.h"
#include "physics/Physics.h"
#include "physics/omega/omega.h"
#include "base/Logger.h"
#include "base/detail/CmdLine.h"

#include "TRint.h"

using namespace std;
using namespace ant::output;
using namespace ant;
using namespace ant::analysis;

int main(int argc, char** argv) {
    SetupLogger(argc, argv);
    TCLAP::CmdLine cmd("Omega Analysis", ' ', "0.1");
    auto input  = cmd.add<TCLAP::MultiArg<string>>("i","input","GoAT input files",true,"string");
    auto output = cmd.add<TCLAP::ValueArg<string>>("o","output","Output file",false,"","string");
    auto batchmode = cmd.add<TCLAP::SwitchArg>("b","batch","Run in batch mode (No ROOT Windows)",false);

    cmd.parse(argc, argv);

    int a=0;
    char** b=nullptr;
    TRint app("omega",&a,b);

    OutputManager om;

    if(output->isSet())
        om.SetNewOutput(output->getValue());

    PhysicsManager pm;

    pm.AddPhysics<OmegaEtaG>(OmegaBase::DataMode::Reconstructed);

    input::GoatReader reader;

    for(auto& file : input->getValue())
            reader.AddInputFile(file);

    reader.Initialize();

    pm.ReadFrom(reader);
    pm.ShowResults();

    app.Run(kTRUE);

    return 0;

}

