#include "AntReader.h"
#include "input/detail/FileManager.h"
#include "TTree.h"
#include <memory>
#include <stdexcept>
#include "detail/Convert.h"
#include "tree/TEvent.h"
#include "data/Event.h"

using namespace std;
using namespace ant;
using namespace ant::input;

AntReader::AntReader():
    files(new FileManager)
{
}

AntReader::~AntReader() {}

void AntReader::AddInputFile(const string& filename)
{
    files->OpenFile(filename);
}


void AntReader::Initialize()
{
    if( !files->GetObject("treeEvent", tree)) {
        throw Exception("Can't find a TTree named treeEvent in input file(s).");
    }


    // Return values of SetBranchAddress:
    // see https://root.cern.ch/root/html534/TTree.html#TTree:CheckBranchAddressType

    if(TTree::kMatch == tree->SetBranchAddress("Event", &buffer)) {
        throw Exception("Can't find a matching branch named \"Event\" in TTree \"treeEvent\"");
    }

}


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
