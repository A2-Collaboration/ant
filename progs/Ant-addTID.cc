#include "simulation/mc/utils/PlutoTID.h"
#include "base/Logger.h"

#include <iostream>

using namespace std;

int main (int argc, char** argv)
{
    SetupLogger();

    if(argc == 2) {
        ant::simulation::mc::utils::PlutoTID::AddTID(std::string(argv[1]));
        return EXIT_SUCCESS ;
    }

    cout << "Add TID Tree to given pluto output file. Will be passed through a2geant (if patched version is used)." << endl;
    cout << "Usage: Ant-addTID <filename>" << endl;

    return EXIT_FAILURE;
}
