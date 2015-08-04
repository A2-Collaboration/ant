#include "OutputManager.h"

#include <iostream>
#include <string>
#include "base/Logger.h"
#include "base/std_ext.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
#include "TFile.h"
#include "TDirectory.h"
#pragma GCC diagnostic pop

using namespace std;
using namespace ant;
using namespace ant::output;



OutputManager::OutputManager() : files()
{
    current_dir = gDirectory;
}

void OutputManager::SetNewOutput(const string &filename)
{

    try {
        auto f = std_ext::make_unique<WrapTFile>(filename,WrapTFile::mode_t::recreate,true);
        current_dir = **f;
        files.emplace_back( std::move(f) );
        VLOG(5) << "Current root output directory is " << current_dir->GetPath();
    } catch (...) {
        cerr << "Can't open output file " << filename << endl;
        current_dir = gDirectory;
    }

}
