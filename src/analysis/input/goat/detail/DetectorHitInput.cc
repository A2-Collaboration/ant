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


    TTree* detectorhit = input_files.GetTree("detectorHits");

    if(detectorhit==nullptr)
        return false;

    NaI.SetupBranches(detectorhit);

//    detectorhit->SetBranchAddress("nNaIHits", &nNaIHits);
//    detectorhit->SetBranchAddress("NaIHits", NaIHits);
//    detectorhit->SetBranchAddress("NaICluster", NaICluster);
//    detectorhit->SetBranchAddress("NaIEnergy", NaIEnergy);
//    detectorhit->SetBranchAddress("NaITime", NaITime);

//    detectorhit->SetBranchAddress("nPIDHits", &nPIDHits);
//    detectorhit->SetBranchAddress("PIDHits", PIDHits);
//    detectorhit->SetBranchAddress("PIDEnergy", PIDEnergy);
//    detectorhit->SetBranchAddress("PIDTime", PIDTime);

//    detectorhit->SetBranchAddress("nMWPCHits", &nMWPCHits);
//    detectorhit->SetBranchAddress("MWPCHits", MWPCHits);

//    detectorhit->SetBranchAddress("nBaF2Hits", &nBaF2Hits);
//    detectorhit->SetBranchAddress("BaF2Hits", BaF2Hits);
//    detectorhit->SetBranchAddress("BaF2Cluster", BaF2Cluster);
//    detectorhit->SetBranchAddress("BaF2Energy", BaF2Energy);
//    detectorhit->SetBranchAddress("BaF2Time", BaF2Time);

//    detectorhit->SetBranchAddress("nVetoHits", &nVetoHits);
//    detectorhit->SetBranchAddress("VetoHits", VetoHits);
//    detectorhit->SetBranchAddress("VetoEnergy", VetoEnergy);
//    detectorhit->SetBranchAddress("VetoTime", VetoTime);

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
