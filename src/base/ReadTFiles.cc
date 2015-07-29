#include "ReadTFiles.h"

#include "base/Logger.h"

using namespace std;
using namespace ant;

ReadTFiles::ReadTFiles() : files()
{

}

ReadTFiles::~ReadTFiles()
{
    CloseAll();
}

bool ReadTFiles::OpenFile(const std::string &filename)
{
    TFile* f = new TFile(filename.c_str(), "READ");

    if(f && !f->IsZombie()) {
        files.emplace_back(f);
        LOG(INFO) << "Opened TFile " << filename << " as input";
        return true;
    }
    LOG(ERROR) << "Could not open TFile " << filename;
    return false;
}

void ReadTFiles::CloseAll()
{
    for(auto& f: files) {
        f->Close();
        LOG(INFO) << "Closed TFile " << f->GetName();
    }

    files.clear();
}

