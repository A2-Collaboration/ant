#pragma once

#include <memory>
#include <string>
#include <list>

//ROOT
#include "TFile.h"
#include "TDirectory.h"
#include "TError.h"
#include "TList.h"
#include "TCollection.h"
#include "TKey.h"

//class TList;

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

    bool isOpen() const;
    bool isZombie() const;


public:

    WrapTFile(const std::string& filename, mode_t access_mode = mode_t::recreate, bool change_gDirectory = false );
    ///@todo is this required?? remove if possible
    TFile* operator* () { return file.get(); }

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

    /**
     * @brief GetListOf Generate a list with all Object in File of provided Type
     * @return list with pointers to object;
     */
    template<class T>
    std::list<T*> GetListOf() const
    {
        std::list<T*> theList;

        TList* keys = GetListOfKeys();

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
        return theList;
    }

    template <typename T>
    void GetObject(const std::string& name, T*& obj) {
        file->GetObject(name.c_str(), obj);
    }


    TList* GetListOfKeys() const
    {
        return file->GetListOfKeys();
    }

    ~WrapTFile();
    WrapTFile(const WrapTFile&) = delete;
    WrapTFile& operator= (const WrapTFile&) = delete;
};

}
