#include "FileManager.h"

#include "TFile.h"

using namespace ant;
using namespace ant::input;

FileManager::FileManager()
{

}

FileManager::~FileManager()
{

}

bool FileManager::OpenFile(const std::string &filename)
{
    TFile* f = new TFile(filename.c_str(), "READ");

    if(f && f->IsOpen()) {
        files.emplace_back(f);
        return true;
    }

    return false;
}

void FileManager::CloseAll()
{
    for(auto& f: files) {
        f->Close();
    }

    files.clear();
}

