#include "WrapTFile.h"

//ant
#include "std_ext.h"
#include "Logger.h"

//ROOT
#include "TFile.h"
#include "TDirectory.h"
#include "TError.h"

//std
#include <stdexcept>
#include <string>

using namespace std;
using namespace ant;

WrapTFile::WrapTFile(const std::string& filename, mode_t access_mode, bool change_gDirectory):
    mode(access_mode),
    changeDirectory(change_gDirectory)
{

    string root_mode;

    switch (mode)
    {
    case mode_t::recreate:
        root_mode = "RECREATE";
        break;
        /*     TODO
    case mode_t::create:
        root_mode = "CREATE";
        break;
    case mode_t::update:
        root_mode = "UPDATE";
        break;
    case mode_t::read:
        root_mode = "READ";
        break;
        */
    default:
        root_mode = "RECREATE";
        break;
    }

    if ( !changeDirectory )
    {
        const auto prev_gErrorIgnoreLevel = gErrorIgnoreLevel;
        const auto prev_Directory         = gDirectory;
        gErrorIgnoreLevel = kError + 1;
        file = std_ext::make_unique<TFile>(filename.c_str(), root_mode.c_str());
        gDirectory = prev_Directory;
        gErrorIgnoreLevel = prev_gErrorIgnoreLevel;
    }
    else
        file = std_ext::make_unique<TFile>(filename.c_str(), root_mode.c_str());

    if(!IsOpen())
        throw std::runtime_error(string("Could not open TFile for writing at ")+filename);
    VLOG(5) << "Opened output file " << filename;
}

bool WrapTFile::IsOpen() const
{
    return file->IsOpen();
}

void WrapTFile::cd()
{
    if(!IsOpen())
        return;
    file->cd();
}

WrapTFile::~WrapTFile()
{
    if(!IsOpen())
        return;

    VLOG(5) << "Syncing output file " << file->GetName();
    file->Write();
    file->Close();
    LOG(INFO) << "Closed output file " << file->GetName();

}
