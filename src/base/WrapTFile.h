#pragma once

#include <memory>
#include <string>

class TFile;

namespace ant {

class WrapTFile {

public:
    enum class mode_t
    {
        create,     //  create new file and open for writing, if it doesn't exist
        recreate,   //  will overwrite any existing file
        update,     //  open existing file for writing, if it doesn't exist create it
        read        //  open file read only
    };

protected:
    std::unique_ptr<TFile> file;
    mode_t mode;
    bool changeDirectory;


public:

    WrapTFile(const std::string& filename, mode_t access_mode = mode_t::recreate, bool change_gDirectory = true );

    ///@todo is this required?? remove if possible
    TFile* operator* () { return file.get(); }

    bool IsOpen() const;
    void cd();

    ~WrapTFile();
    WrapTFile(const WrapTFile&) = delete;
    WrapTFile& operator= (const WrapTFile&) = delete;
};

}
