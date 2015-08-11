#pragma once

#include <memory>
#include <string>
#include <list>

#include "base/std_ext.h"
//ROOT
#include "TFile.h"
#include "TDirectory.h"
#include "TError.h"
#include "TList.h"
#include "TCollection.h"
#include "TKey.h"

class TH1D;
class TH2D;
class TH3D;


namespace ant {



class WrapTFile {
protected:
    std::list<std::unique_ptr<TFile>> files;

    std::unique_ptr<TFile> openFile(const std::string& filename, const std::string mode);

    WrapTFile();

public:

    /**
     * @brief GetListOf Generate a list with all Object in File of provided Type
     * @return list with pointers to object;
     */
    template<class T>
    std::list<T*> GetListOf() const
    {
        std::list<T*> theList;

        for(auto& file : files) {

            TList* keys = file->GetListOfKeys();

            if(!keys)
                return {};

            T* objectPtr = nullptr;
            TKey* key = nullptr;
            TIter nextk(keys);

            while((key = (TKey*)nextk()))
            {
                objectPtr = dynamic_cast<T*>(key->ReadObj());
                if ( !objectPtr )
                    continue;
                theList.push_back(objectPtr);
            }
        }
        return theList;
    }

    template <typename T>
    bool GetObject(const std::string& name, T*& ptr) {
        ptr = nullptr;

        std::string filename;

        for(auto& file : files) {
            T* ptr_ = nullptr;
            file->GetObject(name.c_str(), ptr_);
            if(ptr == nullptr && ptr_ != nullptr) {
                ptr = ptr_;
                filename = file->GetName();
            }
            else if(ptr != nullptr && ptr_ != nullptr) {
                throw std::runtime_error(
                            ant::std_ext::formatter()
                            << "Found object " << name << " in "
                            << filename << " and in " << file->GetName());
            }
        }
        return ptr != nullptr;
    }

    std::shared_ptr<TH1D> GetSharedTH1(const std::string& name);
    std::shared_ptr<TH2D> GetSharedTH2(const std::string& name);
    std::shared_ptr<TH3D> GetSharedTH3(const std::string& name);

    template<typename T>
    std::shared_ptr<T>    GetSharedClone(const std::string& name) {
        T* ptr = nullptr;
        GetObject(name.c_str(), ptr);
        if(ptr) {
            return std::shared_ptr<T>(dynamic_cast<T*>(ptr->Clone()));
        }

        return nullptr;
    }

    virtual ~WrapTFile();
    WrapTFile(const WrapTFile&) = delete;
    WrapTFile& operator= (const WrapTFile&) = delete;
};

class WrapTFileOutput: public WrapTFile {
public:
    enum class mode_t
    {
        create,     //  create new file and open for writing, if it doesn't exist
        recreate,   //  will overwrite any existing file
        update      //  open existing file for writing, if it doesn't exist create it
    };

    WrapTFileOutput(const std::string& filename, mode_t Mode=mode_t::recreate, bool changeDirectory=false);
    virtual ~WrapTFileOutput();

    void cd();

    template<class T, typename... Args>
    T* CreateInside(Args&&... args)
    {
        const auto prev_Directory = gDirectory;
        cd();
        T* object = new T(std::forward<Args>(args)...);
        gDirectory = prev_Directory;

        return object;
    }

    WrapTFileOutput(const WrapTFileOutput&) = delete;
    WrapTFileOutput& operator= (const WrapTFileOutput&) = delete;
};

class WrapTFileInput: public WrapTFile {
public:
    WrapTFileInput();
    WrapTFileInput(const std::string& filename);
    virtual ~WrapTFileInput();

    void OpenFile(const std::string& filename);

    std::size_t NumberOfFiles() const;

    WrapTFileInput(const WrapTFileInput&) = delete;
    WrapTFileInput& operator= (const WrapTFileInput&) = delete;
};

}
