#include "input/goat/GoatReader.h"
#include "OutputManager.h"
#include "physics/Physics.h"
#include "physics/omega/omega.h"

using namespace std;
using namespace ant::output;
using namespace ant;

int main(int argc, char** argv) {

    OutputManager om;

    om.SetNewOutput("omega.root");

    PhysicsManager pm;

    pm.AddPhysics<ant::analysis::Omega3>("Omega3_rec", false);

    input::GoatReader reader;

    for(int i=1; i<argc;++i)
        reader.AddInputFile(argv[i]);

    reader.Initialize();

    pm.ReadFrom(reader);

    return 0;

}

