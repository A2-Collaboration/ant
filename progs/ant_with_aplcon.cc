#include "analysis/physics/test/TestAPLCON.h"
#include "analysis/OutputManager.h"
#include "analysis/input/goat/GoatReader.h"
#include "analysis/physics/Physics.h"

using namespace std;
using namespace ant;
using namespace ant::output;

int main(int argc, char *argv[])
{

    OutputManager om;

    om.SetNewOutput(string(argv[0])+"_output.root");

    PhysicsManager pm;

    pm.AddPhysics<analysis::TestAPLCON>();

    input::GoatReader reader;

    for(int i=1; i<argc;++i)
        reader.AddInputFile(argv[i]);

    reader.Initialize();

    pm.ReadFrom(reader);

    return 0;

}
