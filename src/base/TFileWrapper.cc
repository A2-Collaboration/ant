#include "TFileWrapper.h"
#include "TFile.h"
#include "std_ext.h"
#include "Logger.h"

using namespace std;
using namespace ant;

TFileWrapper::TFileWrapper(const string &filename)
{
    file = std_ext::make_unique<TFile>(filename.c_str(), "RECREATE");

    if(file && file->IsOpen()) {
        VLOG(5) << "Opened output file " << filename;
    } else
        throw false;
}

bool TFileWrapper::isOpen() const
{
    return (file && file->IsOpen());
}

void TFileWrapper::cd()
{
    if(isOpen()) {
        file->cd();
    }
}

TFileWrapper::~TFileWrapper()
{
    if(file) {
        if(file->IsOpen()) {
            VLOG(5) << "Syncing output file " << file->GetName();
            file->Write();
            file->Close();
            VLOG(5) << "Closed output file " << file->GetName();
        }
    }
}
