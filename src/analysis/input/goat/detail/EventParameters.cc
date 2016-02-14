#include "EventParameters.h"

#include "TTree.h"

using namespace ant::analysis::input;


EventParameters::EventParameters()
{

}

EventParameters::~EventParameters()
{

}

bool EventParameters::SetupBranches(TreeRequestManager&& input_files) {

    TTree* tree = input_files.GetTree("eventParameters");

    if(tree==nullptr)
        return false;

    tree->SetBranchAddress("eventNumber",    &EventNumber);
    tree->SetBranchAddress("nReconstructed", &nReconstructed);

    return true;
}

void EventParameters::GetEntry()
{
}
