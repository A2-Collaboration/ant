#include "DetectorHitInput.h"

#include "base/std_ext/string.h"

#include "TTree.h"

#include <iostream>

using namespace std;
using namespace ant::analysis::input;


DetectorHitInput::DetectorHitInput() :
    NaI("NaI"),
    PID("PID"),
    MWPC("MWPC"),
    BaF2("BaF2"),
    Veto("Veto")
{

}

DetectorHitInput::~DetectorHitInput()
{

}

void DetectorHitInput::GetEntry()
{

}

bool DetectorHitInput::SetupBranches(TreeRequestManager&& input_files) {


    TTree* tree = input_files.GetTree("detectorHits");

    if(tree==nullptr)
        return false;

    NaI.SetupBranches(tree);
    PID.SetupBranches(tree);
    BaF2.SetupBranches(tree);
    Veto.SetupBranches(tree);

    return true;

}

DetectorHitInput::Item_t::Item_t(const string& name) :
    Name(name)
{
}

void DetectorHitInput::Item_t::SetupBranches(TTree* tree)
{
    auto branches = tree->GetListOfBranches();
    for(Int_t i=0;i<branches->GetEntries();i++) {
        auto branch = dynamic_cast<TBranch*>((*branches)[i]);
        const string branchname = branch->GetName();
        if(std_ext::string_starts_with(branchname, "n"+Name)) {
            branch->SetAddress(&nHits);
        }
        else if(std_ext::string_starts_with(branchname, Name)) {
            if(std_ext::string_ends_with(branchname, "Hits"))
                branch->SetAddress(Hits);
            if(std_ext::string_ends_with(branchname, "Cluster"))
                branch->SetAddress(Cluster);
            if(std_ext::string_ends_with(branchname, "Energy"))
                branch->SetAddress(Energy);
            if(std_ext::string_ends_with(branchname, "Time"))
                branch->SetAddress(Time);
        }
    }
}
