#include "UnpackerReader.h"

#include "TEvent.h"
#include "TDetectorRead.h"
#include "THeaderInfo.h"
#include "TUnpackerMessage.h"
#include "TSlowControl.h"


#include "base/ReadTFiles.h"
#include "base/Logger.h"

#include <algorithm>


using namespace std;
using namespace ant;
using namespace ant::tree;




UnpackerReader::UnpackerReader(const std::shared_ptr<ReadTFiles>& rootfiles) :
    file(rootfiles), // remember the shared_ptr to make sure it lives as long as this class
    treerecords(),
    isopen(false)
{}

UnpackerReader::~UnpackerReader() {}

bool UnpackerReader::OpenInput() {
    // the order here determines the order of NextItem()
    // if the TID is equal
    if(!SetupBranch("HeaderInfo", HeaderInfo))
        return false;
    if(!SetupBranch("UnpackerMessage", UnpackerMessage))
        return false;
    if(!SetupBranch("SlowControl", SlowControl))
        return false;
    if(!SetupBranch("DetectorRead", DetectorRead))
        return false;
    if(!SetupBranch("Event", Event))
        return false;

    isopen = true;
    return true;
}

bool UnpackerReader::GetUniqueHeaderInfo(THeaderInfo& headerInfo) {
    for(const treerecord_t& treerecord : treerecords) {
        const TDataRecord& rec = treerecord.GetRecord();
        if(rec.IsA() == THeaderInfo::Class()) {
            TTree* tree = treerecord.Tree;
            if(tree->GetEntries() != 1)
                return false;
            tree->GetEntry(0);
            headerInfo = *HeaderInfo;
            return true;
        }
    }
    return false;
}

std::shared_ptr<TDataRecord> UnpackerReader::NextItem() noexcept {

    if(treerecords.empty())
        return nullptr;

    auto compare = [] (const treerecord_t& a, const treerecord_t& b) {
        return a.GetRecord().ID < b.GetRecord().ID;
    };

    auto it_treerecord = min_element(treerecords.begin(), treerecords.end(), compare);

    auto record = shared_ptr<TDataRecord>(*(it_treerecord->Record));
    if(!it_treerecord->GetNext()) {
        // fetched the last item from this treerecord entry,
        // so remove it from the available treerecords
        treerecords.erase(it_treerecord);
    }
    return record;
}
