#include "WrapTFile.h"

#include "std_ext/memory.h"
#include "std_ext/misc.h"
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
    auto errorhandler = GetErrorHandler();
    std_ext::execute_on_destroy restore_error_handler([errorhandler] () {
        SetErrorHandler(errorhandler);
    });
    SetErrorHandler([] (
                    int, Bool_t, const char *location, const char *msg) {
        throw std::runtime_error(std_ext::formatter() << "Could not open TFile: "
                                 << location << ": " << msg);
    });

    auto file = std_ext::make_unique<TFile>(filename.c_str(), mode.c_str());

    if(!file->IsOpen() || file->IsZombie())
    {
        throw std::runtime_error("Could not open TFile at "+filename);
    }

    return file;
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
