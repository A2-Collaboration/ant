#pragma once

#include "TDataRecord.h" // for TID
#include "unpacker/Unpacker.h" // for Unpacker::Reader interface
#include "base/ReadTFiles.h"

#include "TTree.h"

#include <vector>
#include <memory>
#include <string>

namespace ant {

class TEvent;
class TDetectorRead;
class THeaderInfo;
class TUnpackerMessage;
class TSlowControl;

namespace tree {

class UnpackerReader : public Unpacker::Reader {

    std::shared_ptr<ReadTFiles> file;

    TEvent* Event;
    TDetectorRead* DetectorRead;
    THeaderInfo* HeaderInfo;
    TUnpackerMessage* UnpackerMessage;
    TSlowControl* SlowControl;

    struct treerecord_t {
        TDataRecord** Record;
        Long64_t CurrEntry;
        TTree* Tree;
        const TDataRecord& GetRecord() const {
            return *(*Record);
        }
        bool GetNext() {
            ++CurrEntry;
            if(CurrEntry == Tree->GetEntriesFast())
                return false;
            // setting the ptr to zero gives us the
            // control over the allocated data
            *Record = nullptr;
            Tree->GetEntry(CurrEntry);
            return true;
        }
    };

    using treerecords_t = std::list<treerecord_t>;
    treerecords_t treerecords;
    treerecords_t::iterator it_treerecord;
    bool isopen;

    TID currID;

    TID findMinID() const;

public:

    UnpackerReader(const std::shared_ptr<ReadTFiles>& rootfiles);
    virtual ~UnpackerReader();

    template<typename T>
    bool SetupBranch(const std::string& name, T*& ptr) {
        TTree* tree = nullptr;
        const std::string treename = std::string("tree") + name;
        if(!file->GetObject(treename, tree))
            return false;
        if(tree->GetEntries()==0)
            return true;
        if(!tree->GetListOfBranches()->FindObject(name.c_str()))
            return false;
        tree->SetBranchAddress(name.c_str(), std::addressof(ptr));
        treerecords.emplace_back(treerecord_t{
                    reinterpret_cast<TDataRecord**>(std::addressof(ptr)),
                    0,
                    tree}
                    );
        ptr = new T();
        tree->GetEntry(0);
        return true;
    }

    bool OpenInput();

    bool IsOpen() const {
        return isopen;
    }

    bool GetUniqueHeaderInfo(THeaderInfo& headerInfo);

    virtual std::shared_ptr<TDataRecord> NextItem() noexcept override;



//    void Fill(TDataRecord* rec) {
//        const TClass* isA = rec->IsA();
//        if(isA == TEvent::Class()) {
//            Event = static_cast<TEvent*>(rec);
//            treeEvent->Fill();
//        }
//        else if(isA == TDetectorRead::Class()) {
//            DetectorRead = static_cast<TDetectorRead*>(rec);
//            treeDetectorRead->Fill();
//        }
//        else if(isA == THeaderInfo::Class()) {
//            HeaderInfo = static_cast<THeaderInfo*>(rec);
//            treeHeaderInfo->Fill();
//        }
//    }



};

}}
