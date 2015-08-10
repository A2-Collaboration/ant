#include "DetectorHitInput.h"

#include "TTree.h"

using namespace ant::analysis::input;

DetectorHitInput::DetectorHitInput()
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

    detectorhit->SetBranchAddress("nNaIHits", &nNaIHits);
    detectorhit->SetBranchAddress("NaIHits", NaIHits);
    detectorhit->SetBranchAddress("NaICluster", NaICluster);
    detectorhit->SetBranchAddress("nPIDHits", &nPIDHits);
    detectorhit->SetBranchAddress("PIDHits", PIDHits);
    detectorhit->SetBranchAddress("nMWPCHits", &nMWPCHits);
    detectorhit->SetBranchAddress("MWPCHits", MWPCHits);
    detectorhit->SetBranchAddress("nBaF2Hits", &nBaF2Hits);
    detectorhit->SetBranchAddress("BaF2Hits", BaF2Hits);
    detectorhit->SetBranchAddress("BaF2Cluster", BaF2Cluster);
    detectorhit->SetBranchAddress("nVetoHits", &nVetoHits);
    detectorhit->SetBranchAddress("VetoHits", VetoHits);

    return true;

}
