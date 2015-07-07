#include "TriggerInput.h"

#include "TTree.h"

using namespace ant;
using namespace input;



TriggerInput::TriggerInput()
{
    for(Int_t i=0; i<GTreeTrigger_MAX; i++)
    {
        triggerPattern[i] = 0;
        errorModuleID[i] = 0;
        errorModuleIndex[i] = 0;
        errorCode[i] = 0;
    }
}

TriggerInput::~TriggerInput()
{

}


bool TriggerInput::SetupBranches(TreeRequestManager&& input_files) {

    TTree* trigger = input_files.GetTree("trigger");

    if(trigger==nullptr)
        return false;

    trigger->SetBranchAddress("energySum", 	      &energySum);
    trigger->SetBranchAddress("multiplicity", 	  &multiplicity);
    trigger->SetBranchAddress("nTriggerPattern",  &nTriggerPattern);
    trigger->SetBranchAddress("triggerPattern",   triggerPattern);
    trigger->SetBranchAddress("nErrors", 	      &nErrors);
    trigger->SetBranchAddress("errorModuleID",    errorModuleID);
    trigger->SetBranchAddress("errorModuleIndex", errorModuleIndex);
    trigger->SetBranchAddress("errorCode",        errorCode);
    if(trigger->GetBranch("helicity"))
    {
        trigger->SetBranchAddress("helicity",  &helicity);
        hasHelicity = true;
    }
    if(trigger->GetBranch("mc_evt_id"))
    {
        trigger->SetBranchAddress("mc_evt_id", &MC_evt_id);
        trigger->SetBranchAddress("mc_rnd_id", &MC_rnd_id);
        hasMCID = true;
    }

    return true;
}

void TriggerInput::GetEntry()
{
}
