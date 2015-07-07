#include "OutputManager.h"

#include <iostream>
#include <string>
#include "TFile.h"
#include "TDirectory.h"

using namespace std;
using namespace ant;
using namespace ant::output;



OutputManager::OutputManager()
{
    current_dir = gDirectory;
}

void OutputManager::SetNewOutput(const string &filename)
{

    try {
        auto f = std::unique_ptr<TFileWrapper>( new TFileWrapper(filename));
        current_dir = **f;
        files.emplace_back( std::move(f) );
    } catch (...) {
        cerr << "Can't open output file " << filename << endl;
        current_dir = gDirectory;
    }

}

OutputManager::TFileWrapper::TFileWrapper(const string &filename)
{
    TFile* f = new TFile(filename.c_str(), "RECREATE");
    if(f && f->IsOpen())
        file = f;
    else
        throw false;
}

OutputManager::TFileWrapper::~TFileWrapper()
{
    if(file) {
        cout << "closing file " << file->GetName() << endl;
        if(file->IsOpen()) {
            file->Write();
            file->Close();
        }
        delete file;
    }
}
