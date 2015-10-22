#include "WrapTFile.h"

#include "std_ext/memory.h"
#include "Logger.h"

#include "TSystem.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "Compression.h"

#include <stdexcept>
#include <string>

using namespace std;
using namespace ant;

std::unique_ptr<TFile> WrapTFile::openFile(const string& filename, const string mode)
{
    std::unique_ptr<TFile> file;
    const auto prev_gErrorIgnoreLevel = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kError + 1;

    file = std_ext::make_unique<TFile>(filename.c_str(), mode.c_str());

    gErrorIgnoreLevel = prev_gErrorIgnoreLevel;

    if(!file->IsOpen() || file->IsZombie() )
    {
        throw std::runtime_error(string("Could not open TFile at ")+filename);
    }

    return std::move(file);
}

WrapTFile::WrapTFile()
{
}

WrapTFile::~WrapTFile()
{
    for(auto& file : files) {
        file->Close();
        VLOG(5) << "Closed file " << file->GetName();
    }
}

std::shared_ptr<TH1> WrapTFile::GetSharedHist(const string& name)
{
    TH1* hist = nullptr;
    GetObject(name, hist);
    if(hist)
        hist->SetDirectory(nullptr);
    return std::shared_ptr<TH1>(hist);
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




//============================================================================================

struct SavedDirectory_t {
    TDirectory* dir;
    SavedDirectory_t() : dir(gDirectory) {}
    ~SavedDirectory_t() { if(dir) gDirectory = dir; }
    void pop() { gDirectory = dir; dir = nullptr; }
};

WrapTFileOutput::WrapTFileOutput(const std::string& filename, mode_t access_mode, bool changeDirectory)
{

    string root_mode;

    switch (access_mode)
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
    default:
        root_mode = "RECREATE";
        break;
    }

    // recursively create the directory
    stringstream ss_cmd;
    ss_cmd << "mkdir -p " << gSystem->DirName(filename.c_str());
    auto retval = gSystem->Exec(ss_cmd.str().c_str());
    VLOG(5) << "Executed '" << ss_cmd.str() << "' with code " << retval;


    std::unique_ptr<TFile> file;

    if ( !changeDirectory )
    {
        SavedDirectory_t d;
        file = openFile(filename, root_mode);
        d.pop();
    }
    else
        file = openFile(filename, root_mode);

    //file->SetCompressionAlgorithm(ROOT::kLZMA);
    //file->SetCompressionLevel(5);

    files.emplace_back(move(file));

    VLOG(5) << "Opened file " << filename << " in " << root_mode << "-mode.";
}

WrapTFileOutput::~WrapTFileOutput()
{
    LOG(INFO) << "Writing output file " <<  files.front()->GetName();
    files.front()->Write();
}

void WrapTFileOutput::cd()
{
    files.front()->cd();
}



//============================================================================================



WrapTFileInput::WrapTFileInput()
{
}

WrapTFileInput::WrapTFileInput(const string& filename)
{
    OpenFile(filename);
}

WrapTFileInput::~WrapTFileInput()
{
}

void WrapTFileInput::OpenFile(const string& filename)
{
    SavedDirectory_t d;
    auto file = openFile(filename, "READ");
    VLOG(5) << "Opened file " << filename << " for reading";
    d.pop();

    files.emplace_back(std::move(file));
}

size_t WrapTFileInput::NumberOfFiles() const
{
    return files.size();
}
