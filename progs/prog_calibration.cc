#include "analysis/input/goat/GoatReader.h"
#include "analysis/OutputManager.h"
#include "analysis/physics/Physics.h"
#include "analysis/physics/omega/omega.h"
#include "base/Logger.h"
#include "calibration/BaseCalModule.h"
#include "calibration/TestCalCB.h"
#include <string>

using namespace std;
using namespace ant::output;
using namespace ant;

int main(int argc, char** argv) {
    SetupLogger(argc, argv);

    OutputManager om;

    om.SetNewOutput("calibration.root");

    PhysicsManager pm;

    pm.AddPhysics<ant::calibration::TestCalCB>();

    input::GoatReader reader;

    for(int i=1; i<argc;++i) {
        if(argv[i][0] != '-')
            reader.AddInputFile(argv[i]);
    }

    reader.Initialize();

    pm.ReadFrom(reader);

    return 0;

}

