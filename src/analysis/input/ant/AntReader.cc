#include "AntReader.h"

#include "data/Event.h"
#include "detail/Convert.h"
#include "base/ReadTFiles.h"
#include "base/std_ext.h"
#include "tree/TEvent.h"

#include "TTree.h"

#include <memory>
#include <stdexcept>

using namespace std;
using namespace ant;
using namespace ant::input;

AntReader::AntReader(const std::shared_ptr<ReadTFiles>& rootfiles) :
    files(rootfiles)
{
    if( !files->GetObject("treeEvent", tree)) {
        throw Exception("Can't find a TTree named treeEvent in input file(s).");
    }


    // Return values of SetBranchAddress:
    // see https://root.cern.ch/root/html534/TTree.html#TTree:CheckBranchAddressType

    const auto result = tree->SetBranchAddress("Event", &buffer);

    if(result != TTree::kMatch) {
        throw Exception("Can't find a matching branch named \"Event\" in TTree \"treeEvent\" (" + to_string(result) + ")");
    }


}

AntReader::~AntReader() {}

Long64_t AntReader::GetNEvents() const
{
    return tree->GetEntries();
}


std::shared_ptr<Event> AntReader::ReadNextEvent()
{
    tree->GetEntry(++current);

    auto event = input::Convert(*buffer);

    return event;
}

bool AntReader::hasData() const {
    return current < GetNEvents();
}

long long AntReader::EventsRead() const
{
    return current+1;
}

long long AntReader::TotalEvents() const
{
    return GetNEvents();
}
