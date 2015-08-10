#include "WrapTFile.h"

//ant
#include "std_ext.h"
#include "Logger.h"


//std
#include <stdexcept>
#include <string>

#include "TSystem.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"

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
    case mode_t::create:
        root_mode = "CREATE";
        break;
    case mode_t::update:
        root_mode = "UPDATE";
        break;
    case mode_t::read:
        root_mode = "READ";
        break;
    default:
        root_mode = "RECREATE";
        break;
    }

    if(mode != mode_t::read) {
        // recursively create the directory
        stringstream ss_cmd;
        ss_cmd << "mkdir -p " << gSystem->DirName(filename.c_str());
        VLOG(5) << "Executed '" << ss_cmd.str() << "'' with code "
                << gSystem->Exec(ss_cmd.str().c_str());
    }

    const auto prev_gErrorIgnoreLevel = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kError + 1;

    if ( !changeDirectory )
    {
        const auto prev_Directory         = gDirectory;
        file = std_ext::make_unique<TFile>(filename.c_str(), root_mode.c_str());
        gDirectory = prev_Directory;
    }
    else
        file = std_ext::make_unique<TFile>(filename.c_str(), root_mode.c_str());

    gErrorIgnoreLevel = prev_gErrorIgnoreLevel;


    if(!isOpen() || isZombie() )
    {
        if ( mode == mode_t::read )
            throw runtime_error(std_ext::formatter() << "Could not open " << filename << " for reading.");
        else
            throw std::runtime_error(string("Could not open TFile for writing at ")+filename);
    }

    VLOG(5) << "Opened file " << filename << " in " << root_mode << "-mode.";
}

bool WrapTFile::isOpen() const
{
    return file->IsOpen();
}

bool WrapTFile::isZombie() const
{
    return file->IsZombie();
}

void WrapTFile::cd()
{
    if(!isOpen())
        return;
    file->cd();
}

std::shared_ptr<TH1D> WrapTFile::GetSharedTH1(const string& name)
{
    TH1D* hist = nullptr;
    GetObject(name, hist);
    if(hist)
        hist->SetDirectory(nullptr);
    return std::shared_ptr<TH1D>(hist);
}

std::shared_ptr<TH2D> WrapTFile::GetSharedTH2(const string& name)
{
    TH2D* hist = nullptr;
    GetObject(name, hist);
    if(hist)
        hist->SetDirectory(nullptr);
    return std::shared_ptr<TH2D>(hist);
}

std::shared_ptr<TH3D> WrapTFile::GetSharedTH3(const string& name)
{
    TH3D* hist = nullptr;
    GetObject(name, hist);
    if(hist)
        hist->SetDirectory(nullptr);
    return std::shared_ptr<TH3D>(hist);
}

WrapTFile::~WrapTFile()
{
    if(!isOpen())
        return;

    if (mode == mode_t::read)
    {
        file->Close();
        LOG(INFO) << "Closed file " << file->GetName();
        return;
    }

    VLOG(5) << "Syncing output file " << file->GetName();
    file->Write();
    file->Close();
    LOG(INFO) << "Closed output file " << file->GetName();

}
