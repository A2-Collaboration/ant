#include "ReadTFiles.h"

#include "base/Logger.h"
#include "base/std_ext.h"

#include "TError.h"

using namespace std;
using namespace ant;

ReadTFiles::ReadTFiles() : files() {}

ReadTFiles::~ReadTFiles()
{
    CloseAll();
}

bool ReadTFiles::OpenFile(const std::string& filename, bool silent)
{
    const auto prev_gErrorIgnoreLevel = gErrorIgnoreLevel;
    if(silent)
        gErrorIgnoreLevel = kError+1;
    auto tfile = std_ext::make_unique<TFile>(filename.c_str(), "READ");
    gErrorIgnoreLevel = prev_gErrorIgnoreLevel;

    if(tfile->IsZombie()) {
        if(!silent)
            LOG(ERROR) << "Could not open TFile " << filename;
        return false;
    }

    files.emplace_back(move(tfile));
    return true;
}

void ReadTFiles::CloseAll()
{
    for(auto& f: files) {
        f->Close();
        LOG(INFO) << "Closed TFile " << f->GetName();
    }

    files.clear();
}

