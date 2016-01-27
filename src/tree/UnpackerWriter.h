#pragma once

#include "unpacker/Unpacker.h" // for Unpacker::Reader interface

#include "tree/MemoryPool.h"
#include "tree/TDataRecord.h"

#include "TTree.h"

#include <vector>
#include <memory>
#include <string>

namespace ant {

struct TEvent;
struct TDetectorRead;
struct TUnpackerMessage;
struct TSlowControl;
class WrapTFile;

namespace tree {

class UnpackerWriter {

    std::unique_ptr<WrapTFile> file;

    template<class T>
    struct treerecord_t {
        treerecord_t() :
            Record(nullptr),
            Tree(nullptr)
        {}

        // call this after opening the TFile
        void Init() {
            const std::string classname = T::Class()->GetName();
            size_t Tpos = classname.find_first_of('T');
            if(Tpos == std::string::npos)
                return;
            Tpos++; // get past the T
            const std::string& branchname = classname.substr(Tpos);
            if(branchname.empty())
                return;
            const std::string treename = std::string("tree")+branchname;

            Tree = new TTree(treename.c_str(), treename.c_str());
            Record = new T();
            Tree->Branch(branchname.c_str(), classname.c_str(), &Record);
        }

        bool TryFill(TDataRecord* record) {
            if(Tree==nullptr)
                return false;
            if(record->IsA() != T::Class())
                return false;
            Record = reinterpret_cast<T*>(record);
            Tree->Fill();
            return true;
        }
    private:
        T* Record;
        TTree* Tree;
    };

    treerecord_t<TEvent> Event;
    treerecord_t<TDetectorRead> DetectorRead;
    treerecord_t<TUnpackerMessage> UnpackerMessage;
    treerecord_t<TSlowControl> SlowControl;

public:

    UnpackerWriter(const std::string& outputfile);
    virtual ~UnpackerWriter();

    void Fill(TDataRecord* record) noexcept;

};

}}
