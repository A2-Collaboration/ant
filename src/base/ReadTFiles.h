#pragma once

#include "base/std_ext.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
#include "TFile.h"
#pragma GCC diagnostic pop

#include <list>
#include <string>
#include <memory>
#include <stdexcept>


namespace ant {

class ReadTFiles {
protected:
    using file_list_t = std::list<std::unique_ptr<TFile>>;
    file_list_t files;

public:

    ReadTFiles();
    virtual ~ReadTFiles();

    bool OpenFile(const std::string& filename);

    template <typename T>
    bool GetObject(const std::string& name, T*& ptr) {
        ptr = nullptr;

        std::string filename;

        for(auto& f : files) {
            T* ptr_ = nullptr;
            f->GetObject(name.c_str(), ptr_);
            if(ptr == nullptr && ptr_ != nullptr) {
                ptr = ptr_;
                filename = f->GetName();
            }
            else if(ptr != nullptr && ptr_ != nullptr) {
                throw std::runtime_error(
                            std_ext::formatter()
                            << "Found object " << name << " in "
                            << filename << " and in " << f->GetName());
            }
        }
        return ptr != nullptr;
    }

    void CloseAll();

};

}

