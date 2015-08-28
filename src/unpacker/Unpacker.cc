#include "Unpacker.h"
#include "UnpackerAcqu.h"
#include "UnpackerA2Geant.h"

#include "base/Logger.h"

#include <algorithm>
#include <iostream>


using namespace std;
using namespace ant;

std::unique_ptr<Unpacker::Reader> Unpacker::Get(const string& filename)
{
    // make a list of available unpackers
    std::list< std::unique_ptr<Module> > modules;
    modules.push_back(std_ext::make_unique<UnpackerAcqu>());
    modules.push_back(std_ext::make_unique<UnpackerA2Geant>());

    // remove the unpacker if it says that it could not open the file
    modules.remove_if([&filename] (const unique_ptr<Module>& m) {
        return !m->OpenFile(filename);
    });

    // check if something reasonable is left
    if(modules.empty()) {
        throw Exception("No suitable unpacker found for file "+filename);
    }
    if(modules.size()>1) {
        throw Exception("More than one unpacker found for file "+filename);
    }

    // hand over the unique ptr
    return std::move(modules.back());
}



