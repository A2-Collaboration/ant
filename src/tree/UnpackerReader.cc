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
    it_treerecord(treerecords.end()), // invalid iterator position
    isopen(false),
    currID()
{}

UnpackerReader::~UnpackerReader() {}

bool UnpackerReader::OpenInput() {
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

    if(treerecords.empty())
        return false;

    // set it_treerecord to the item with the minimum ID
    it_treerecord = treerecords.begin();
    currID = findMinID();
    while(it_treerecord != treerecords.end()) {
        if(currID == it_treerecord->GetRecord().ID)
            break;
        ++it_treerecord;
    }
    if(it_treerecord == treerecords.end())
        return false;

    VLOG(9) << "Start reading at ID=" << currID << " with type="
            << it_treerecord->GetRecord().IsA()->GetName();

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

TID UnpackerReader::findMinID() const {
    TID tid_min = it_treerecord->GetRecord().ID; // start with something existing
    // search all for minimum TID
    for(const treerecord_t& treerecord : treerecords) {
        const TDataRecord& record = treerecord.GetRecord();
        if(record.ID < tid_min)
            tid_min = record.ID;
    }
    return tid_min;
}

std::shared_ptr<TDataRecord> UnpackerReader::NextItem() noexcept {

    if(it_treerecord == treerecords.end())
        return nullptr;

    auto record = shared_ptr<TDataRecord>(*(it_treerecord->Record));
    if(!it_treerecord->GetNext()) {
        it_treerecord = treerecords.erase(it_treerecord);
        if(treerecords.empty())
            return record;
        //advance(it_treerecord, -1);
    }

    // find the global minimum
    const TID& tid_min = findMinID();

    // one left anyway
    if(treerecords.size()==1)
        return record;

    // go to other
    do {
        ++it_treerecord;
        if(it_treerecord==treerecords.end())
            it_treerecord = treerecords.begin();

        if(it_treerecord->GetRecord().ID == tid_min)
            return record;
    }
    while(it_treerecord != treerecords.end());

    // should never reach this point
    return nullptr;

}
