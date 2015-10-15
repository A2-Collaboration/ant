#pragma once

#include <string>

namespace ant {
namespace simulation
{
namespace mc
{
namespace utils
{

struct PlutoTID {

    static const std::string tidtree_name;
    static const std::string geant_tidtree_name;

    /**
     * @brief Add a TID Tree to a pluto generated ROOT file.
     * @param filename File to edit
     *
     * Opens the ROOT file in read/write, looks for a "data" TTree and then adds a TID in a new TTree
     * called "dataTID" for each entry in "data"
     */
    static void AddTID(const std::string& filename);

    static void CopyTIDPlutoGeant(const std::string& pluto_filename, const std::string& geant_filename);
};

}
}
}
}
