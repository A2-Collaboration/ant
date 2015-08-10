#include "TaggerInput.h"

#include "TTree.h"

using namespace ant::analysis::input;

TaggerInput::TaggerInput()
{

}

TaggerInput::~TaggerInput()
{

}

void TaggerInput::GetEntry()
{

}

bool TaggerInput::SetupBranches(TreeRequestManager&& input_files){

    TTree* tagger = input_files.GetTree("tagger");

    if(tagger==nullptr)
        return false;

    tagger->SetBranchAddress("nTagged", 	  &nTagged);
    tagger->SetBranchAddress("taggedChannel",  taggedChannel);
    tagger->SetBranchAddress("taggedTime",     taggedTime);
    if(tagger->GetBranch("taggedEnergy"))
    {
        tagger->SetBranchAddress("taggedEnergy",  taggedEnergy);
        hasEnergy = true;
    }

    return true;
}
