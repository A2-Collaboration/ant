#include "UnpackerReader.h"

#include "TEvent.h"
#include "TUnpackerMessage.h"
#include "TSlowControl.h"


#include "base/WrapTFile.h"
#include "base/Logger.h"

#include <algorithm>


using namespace std;
using namespace ant;
using namespace ant::tree;




UnpackerReader::UnpackerReader(const std::shared_ptr<WrapTFileInput>& rootfiles) :
    files(rootfiles), // remember the shared_ptr to make sure it lives as long as this class
    treerecords(),
    isopen(false)
{}

UnpackerReader::~UnpackerReader() {}

double UnpackerReader::PercentDone() const
{
    return double(entries_read)/double(total_entries);
}

bool UnpackerReader::OpenInput() {
    // the order here determines the order of NextItem()
    // if the TID is equal
    if(!SetupBranch("UnpackerMessage", UnpackerMessage))
        return false;
    if(!SetupBranch("SlowControl", SlowControl))
        return false;
    if(!SetupBranch("Event", Event))
        return false;

    isopen = true;
    return true;
}

std::unique_ptr<TDataRecord> UnpackerReader::NextItem() noexcept {

    if(treerecords.empty())
        return nullptr;

    auto compare = [] (const treerecord_t& a, const treerecord_t& b) {
        return a.GetRecord().ID < b.GetRecord().ID;
    };

    auto it_treerecord = min_element(treerecords.begin(), treerecords.end(), compare);

    auto record = unique_ptr<TDataRecord>(*(it_treerecord->Record));
    if(!it_treerecord->GetNext()) {
        // fetched the last item from this treerecord entry,
        // so remove it from the available treerecords
        treerecords.erase(it_treerecord);
    }
    ++entries_read;
    return record;
}
