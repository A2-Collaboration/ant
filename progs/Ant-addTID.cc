#include "simulation/mc/utils/PlutoTID.h"
#include "base/Logger.h"

int main (int argc, char** argv)
{
    SetupLogger();

    if(argc == 2) {
        ant::simulation::mc::utils::PlutoTID::AddTID(std::string(argv[1]));
        return EXIT_SUCCESS ;
    }

    return EXIT_FAILURE;
}
