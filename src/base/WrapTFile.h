#pragma once

#include "base/std_ext/string.h"
#include "base/std_ext/misc.h"

//ROOT
#include "TFile.h"
#include "TDirectory.h"
#include "TError.h"
#include "TList.h"
#include "TCollection.h"
#include "TKey.h"

#include <memory>
#include <string>
#include <list>
#include <set>
#include <stdexcept>
#include <functional>

class TH1;
class TH1D;
class TH2D;
class TH3D;


namespace ant {



class WrapTFile {
protected:
    std::list<std::unique_ptr<TFile>> files;
    std::unique_ptr<TFile> openFile(const std::string& filename, const std::string mode);
    WrapTFile(); // cannot be directly constructed, use WrapTFileInput or WrapTFileOutput

    static bool hasROOTmagic(const std::string& filename);
public:

    /**
     * @brief Traverse applies given function to each leaf in the tree of TDirectory
     * @param func use TKey to inspect anything of that object
     * @note A TKey may exist more than once under the same name but with different cycle number
     */
    void Traverse(std::function<void(TKey*)> func);

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

            // first create a unique list of names,
            // this prevents object with different cycles
            std::set<std::string> names;
            TIter nextk(keys);
            while(auto key = dynamic_cast<TKey*>(nextk()))
            {
                names.emplace(key->GetName());
            }

            for(auto& name : names) {
                T* objectPtr = nullptr;
                file->GetObject(name.c_str(), objectPtr);
                if(objectPtr != nullptr)
                    theList.push_back(objectPtr);
            }
        }
        return theList;
    }

    /**
     * @brief Find TObject of type T and given name in the file(s).
     * @param name Name of the object to get
     * @param ptr  Reference where to store the pointer to the found object
     * @return true if object was found, false otherwise
     */
    template <typename T>
    bool GetObject(const std::string& name, T*& ptr) const {
        ptr = nullptr;

        std::string filename;

        for(auto& file : files) {
            T* ptr_ = nullptr;
            file->GetObject(name.c_str(), ptr_);
            if(ptr_ == nullptr) {
                // GetObject did not work, try GetKey fallback
                // useful for objects which have / in their names
                TKey* key = file->GetKey(name.c_str());
                if(key != nullptr) {
                    ptr_ = dynamic_cast<T*>(key->ReadObj());
                }
            }

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

    template<typename Hist>
    std::shared_ptr<Hist>  GetSharedHist(const std::string& name) const {
        Hist* hist = nullptr;
        GetObject(name, hist);
        if(hist)
            hist->SetDirectory(nullptr);
        return std::shared_ptr<Hist>(hist);
    }


    template<typename T>
    std::shared_ptr<T>  GetSharedClone(const std::string& name) const {
        T* ptr = nullptr;
        if(GetObject(name.c_str(), ptr)) {
            return std::shared_ptr<T>(dynamic_cast<T*>(ptr->Clone()));
        }
        return nullptr;
    }

    template<typename T>
    bool GetObjectClone(const std::string& name, T& object) const {
        T* ptr = nullptr;
        if(GetObject(name.c_str(), ptr)) {
            object = std::move(*ptr);
            return true;
        }
        return false;
    }

    struct Exception : std::runtime_error {
        using std::runtime_error::runtime_error;
    };

    struct ENotARootFile : public Exception {
        using Exception::Exception;
    };


    struct ENotReadable : public Exception {
        using Exception::Exception;
    };

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
        std_ext::execute_on_destroy restoreDir(
                    [prev_Directory] () {
            gDirectory = prev_Directory;
        });

        cd();
        T* object = new T(std::forward<Args>(args)...);

        return object;
    }

    template <typename T>
    int WriteObject(T* obj, const std::string& name) {
        return files.front()->WriteObject(obj, name.c_str());
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

    std::string FileNames() const;

    WrapTFileInput(const WrapTFileInput&) = delete;
    WrapTFileInput& operator= (const WrapTFileInput&) = delete;


};

}
