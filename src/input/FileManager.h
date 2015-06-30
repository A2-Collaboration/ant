#ifndef FILEMANAGER_H
#define FILEMANAGER_H

#include <list>
#include <string>
#include <memory>

#include "TFile.h"

namespace ant {
namespace input {

class FileManager {
protected:
    using file_list_t = std::list<std::unique_ptr<TFile>>;

    file_list_t files;

public:

    FileManager();
    virtual ~FileManager();

    virtual bool OpenFile(const std::string& filename);

    template <typename T>
    bool GetObject(const std::string& name, T*& ptr) {
        ptr = nullptr;

        for(auto& f : files) {
            f->GetObject(name.c_str(), ptr);
            if(ptr) break;
        }

        return ptr != nullptr;
    }

    virtual void CloseAll();


};

}
}

#endif
