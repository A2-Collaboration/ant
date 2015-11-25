#pragma once

#include "unpacker/Unpacker.h" // for Unpacker::Reader interface
#include "base/WrapTFile.h"
#include "base/std_ext/convert.h"
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

/**
 * @brief The UnpackerReader class reads synchronized data from
 */
class UnpackerReader : public Unpacker::Reader {

    std::shared_ptr<WrapTFileInput> files;

    TEvent* Event;
    TDetectorRead* DetectorRead;
    THeaderInfo* HeaderInfo;
    TUnpackerMessage* UnpackerMessage;
    TSlowControl* SlowControl;

    unsigned long long entries_read  =0;
    unsigned long long total_entries =0;

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
    bool isopen;

    template<typename T>
    bool SetupBranch(const std::string& name, T*& ptr) {
        TTree* tree = nullptr;
        ptr = nullptr;
        const std::string treename = std::string("tree") + name;
        if(!files->GetObject(treename, tree))
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
        total_entries += std_ext::clipNegative<decltype(total_entries)>(tree->GetEntries());
        return true;
    }

public:

    UnpackerReader(const std::shared_ptr<WrapTFileInput>& rootfiles);
    virtual ~UnpackerReader();

    double PercentDone() const override;



    bool OpenInput();
    bool IsOpen() const { return isopen; }

    bool GetUniqueHeaderInfo(THeaderInfo& headerInfo);
    virtual std::unique_ptr<TDataRecord> NextItem() noexcept override;

};

}}
