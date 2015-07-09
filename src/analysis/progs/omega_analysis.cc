#include "input/goat/GoatReader.h"
#include "OutputManager.h"
#include "physics/Physics.h"
#include "physics/omega/omega.h"
#include "base/Logger.h"

using namespace std;
using namespace ant::output;
using namespace ant;
using namespace ant::analysis;

int main(int argc, char** argv) {
    SetupLogger(argc, argv);

    OutputManager om;

    om.SetNewOutput("omega.root");

    PhysicsManager pm;

    pm.AddPhysics<OmegaEtaG>(OmegaBase::DataMode::Reconstructed);

    input::GoatReader reader;

    for(int i=1; i<argc;++i) {
        if(argv[i][0] != '-')
            reader.AddInputFile(argv[i]);
    }

    reader.Initialize();

    pm.ReadFrom(reader);

    return 0;

}

