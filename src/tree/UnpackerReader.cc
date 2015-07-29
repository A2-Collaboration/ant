#include "UnpackerReader.h"

#include "TEvent.h"
#include "TDetectorRead.h"
#include "THeaderInfo.h"
#include "TUnpackerMessage.h"
#include "TSlowControl.h"

#include "base/ReadTFiles.h"
#include "base/Logger.h"

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
    if(!SetupBranch("Event", Event))
        return false;
    if(!SetupBranch("DetectorRead", DetectorRead))
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

std::shared_ptr<TDataRecord> UnpackerReader::NextItem() noexcept {
    // the TID determines, what item is next,
    // if there are multiple, the current iterator decides


    return nullptr;

}
