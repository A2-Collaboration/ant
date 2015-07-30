#include "WrapTFile.h"
#include "TFile.h"
#include "std_ext.h"
#include "Logger.h"

#include <stdexcept>

using namespace std;
using namespace ant;

WrapTFile::WrapTFile(const string &filename)
{
    file = std_ext::make_unique<TFile>(filename.c_str(), "RECREATE");
    if(!isOpen())
        throw std::runtime_error(string("Could not open TFile for writing at ")+filename);
    VLOG(5) << "Opened output file " << filename;
}

bool WrapTFile::isOpen() const
{
    return file->IsOpen();
}

void WrapTFile::cd()
{
    if(!isOpen())
        return;
    file->cd();
}

WrapTFile::~WrapTFile()
{
    if(!isOpen())
        return;

    VLOG(5) << "Syncing output file " << file->GetName();
    file->Write();
    file->Close();
    VLOG(5) << "Closed output file " << file->GetName();

}
