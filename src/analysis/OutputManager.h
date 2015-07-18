#ifndef OUTPUTMANAGER_H
#define OUTPUTMANAGER_H

#include <list>
#include <memory>
#include <string>

class TFile;
class TDirectory;

namespace ant {
namespace output {

class OutputManager {
protected:

    class TFileWrapper {
    protected:
        TFile* file = nullptr;

    public:
        TFileWrapper(const std::string& filename);

        TFile* operator* () { return file; }

        ~TFileWrapper();
        TFileWrapper(const TFileWrapper&) = delete;
        TFileWrapper& operator= (const TFileWrapper&) = delete;
    };

    using file_list_t = std::list< std::unique_ptr< TFileWrapper > >;
    file_list_t files;
    TDirectory* current_dir = nullptr;

public:
    OutputManager();

    void SetNewOutput(const std::string& filename);

    TDirectory* CurrentDirectory() { return current_dir; }

    OutputManager(const OutputManager&) = delete;
    OutputManager& operator= (const OutputManager&) = delete;
};

}
}

#endif
