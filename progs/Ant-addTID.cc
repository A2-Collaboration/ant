/**
  * @file Ant-addTID.cc
  * @brief Program to add TID tree to a pluto ROOT file or a geant ROOT file.
  *
  * Adds a TTree to the given pluto ROOT file, containing a branch of TID with one ID for each pluto event (entry in the "data" TTree).
  *
  * The "copy-to-geant" mode copies a TID TTree from a pluto ROOT file to a geant ROOT file.
  * Useful if the TID tree was not added after the pluto file was created and geant has aleady simulated the pluto data.
  * This should be used with caution and if you are very sure that the pluto and geant file belong together.
  */

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
    } else if(argc == 4 && string(argv[1]) == "--copy-to-geant") {
        ant::simulation::mc::utils::PlutoTID::CopyTIDPlutoGeant(argv[2], argv[3]);
        return EXIT_SUCCESS;
    }

    cout << "Add TID Tree to given pluto output file. Will be passed through a2geant (if patched version is used)." << endl;
    cout << "Usage: Ant-addTID <filename>" << endl;
    cout << "Usage: Ant-addTID --copy-to-geant <pluto-file-with-tid> <geant-file-without-tid>" << endl;

    return EXIT_FAILURE;
}
