#include "Physics.h"
#include "OutputManager.h"
#include "input/GoatReader.h"
#include "analysis/DataOverview.h"

using namespace std;
using namespace ant;
using namespace ant::output;

int main(int argc, char *argv[])
{

    OutputManager om;

    om.SetNewOutput("data_overview.root");

    PhysicsManager pm;

    pm.AddPhysics<analysis::DataOverview>();

    input::GoatReader reader;

    for(int i=1; i<argc;++i)
        reader.AddInputFile(argv[i]);

    reader.Initialize();

    pm.ReadFrom(reader);

    return 0;
}

