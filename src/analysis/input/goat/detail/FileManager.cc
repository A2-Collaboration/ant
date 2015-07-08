#include "FileManager.h"

#include "TFile.h"
#include "base/Logger.h"

using namespace ant;
using namespace ant::input;

FileManager::FileManager()
{

}

FileManager::~FileManager()
{
    CloseAll();
}

bool FileManager::OpenFile(const std::string &filename)
{
    TFile* f = new TFile(filename.c_str(), "READ");

    if(f && f->IsOpen()) {
        files.emplace_back(f);
        VLOG(5) << "Opened GoAT File " << filename;
        return true;
    }
    LOG(ERROR) << "Could not open GoAT File " << filename;
    return false;
}

void FileManager::CloseAll()
{
    for(auto& f: files) {
        f->Close();
        VLOG(5) << "Closed GoAT File " << f->GetName();
    }

    files.clear();
}

