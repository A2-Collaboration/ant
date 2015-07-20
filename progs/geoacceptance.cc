#include "analysis/input/goat/GoatReader.h"
#include "analysis/OutputManager.h"
#include "analysis/physics/Physics.h"
#include "analysis/physics/omega/GeoAcceptance.h"

using namespace std;
using namespace ant::output;
using namespace ant;

int main(int argc, char** argv) {

    OutputManager om;

    om.SetNewOutput("geoacceptance.root");

    PhysicsManager pm;

    pm.AddPhysics<analysis::GeoAcceptance>();

    input::GoatReader reader;

    for(int i=1; i<argc;++i)
        reader.AddInputFile(argv[i]);

    reader.Initialize();

    pm.ReadFrom(reader);

    return 0;

}

