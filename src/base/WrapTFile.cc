#include "WrapTFile.h"

#include "std_ext/memory.h"
#include "std_ext/misc.h"
#include "std_ext/system.h"
#include "std_ext/string.h"
#include "Logger.h"

#include "TSystem.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "Compression.h"
#include "TClass.h"

#include <stdexcept>
#include <string>

using namespace std;
using namespace ant;

static bool warningDetected = false; // need that ugly global variable for ROOT's error handler

std::unique_ptr<TFile> WrapTFile::openFile(const string& filename, const string mode)
{
    auto prevErrorHandler = GetErrorHandler();
    std_ext::execute_on_destroy restore_error_handler([prevErrorHandler] () {
        SetErrorHandler(prevErrorHandler);
    });
    warningDetected = false;
    SetErrorHandler([] (int, Bool_t, const char*, const char*) {
        warningDetected = true;
    });

    auto file = std_ext::make_unique<TFile>(filename.c_str(), mode.c_str());

    if(!file->IsOpen() || file->IsZombie()) {
        throw Exception("Could not properly open TFile at "+filename);
    }

    if(warningDetected) {
        throw Exception("Warning(s) detected when opening "+filename);
    }

    return file;
}

WrapTFile::WrapTFile()
{
}

bool WrapTFile::hasROOTmagic(const std::string& filename)
{
    ifstream is(filename.c_str());
    string buffer(4,'\0');
    is.read(addressof(buffer[0]), 4);
    return buffer == "root";
}

void traverse_dir(TDirectory* dir, std::function<void (TKey*)> func) {

    TList* keys = dir->GetListOfKeys();
    if(!keys)
        return;

    TIter nextk(keys);
    TKey* key;
    while((key = (TKey*)nextk()))
    {
        TClass *classPtr = TClass::GetClass(key->GetClassName());
        if(classPtr->InheritsFrom(TDirectory::Class())) {
            auto subdir = dynamic_cast<TDirectory*>(key->ReadObj());
            if(subdir)
                traverse_dir(subdir, func);
            continue;
        }
        func(key);
    }
}

void WrapTFile::Traverse(std::function<void (TKey*)> func) {
    for(auto& file : files) {
        traverse_dir(file.get(), func);
    }
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
    files.front()->Write();
    LOG(INFO) << "Wrote output file " <<  files.front()->GetName()
              << " (" << (double)files.front()->GetBytesWritten()/(1 << 20) << " MB)";
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
    // test some basic things before giving the file to ROOT
    string errmsg;
    if(!std_ext::system::testopen(filename, errmsg))
        throw ENotReadable(filename+" cannot be opened: "+errmsg);
    if(!hasROOTmagic(filename))
        throw ENotARootFile(filename+" is not a ROOT file");

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

string WrapTFileInput::FileNames() const
{
    std_ext::formatter s;
    s << "(";
    for (const auto& f :files) {
        s << f->GetName() << ", ";
    }
    s<< ")";
    return s.str();
}
