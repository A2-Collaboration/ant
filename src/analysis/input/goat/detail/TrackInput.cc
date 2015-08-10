#include "TrackInput.h"
#include "TTree.h"

using namespace ant::analysis::input;

TrackInput::TrackInput()
{
    for(Int_t i=0; i<GTreeTrack_MAX; i++)
    {
        clusterEnergy[i] = 0;
        theta[i] = 0;
        phi[i] = 0;
        time[i] = 0;
        clusterSize[i] = 0;
        centralCrystal[i] = -1;
        centralVeto[i] = -1;
        detectors[i] = 0;
        vetoEnergy[i] = 0;
        MWPC0Energy[i] = 0;
        MWPC1Energy[i] = 0;
        shortEnergy[i] = 0;
        pseudoVertexX[i] = 0;
        pseudoVertexY[i] = 0;
        pseudoVertexZ[i] = 0;
    }
}

TrackInput::~TrackInput()
{

}

bool TrackInput::SetupBranches(TreeRequestManager &&input_files) {

    TTree* tracks = input_files.GetTree("tracks");

    if(tracks == nullptr)
        return false;

    tracks->SetBranchAddress("nTracks",&nTracks);
    tracks->SetBranchAddress("clusterEnergy",  clusterEnergy);
    tracks->SetBranchAddress("theta", theta);
    tracks->SetBranchAddress("phi",  phi);
    tracks->SetBranchAddress("time", time);
    tracks->SetBranchAddress("clusterSize", clusterSize);
    tracks->SetBranchAddress("centralCrystal", centralCrystal);
    tracks->SetBranchAddress("centralVeto", centralVeto);
    tracks->SetBranchAddress("detectors", detectors);
    tracks->SetBranchAddress("vetoEnergy", vetoEnergy);
    tracks->SetBranchAddress("MWPC0Energy", MWPC0Energy);
    tracks->SetBranchAddress("MWPC1Energy", MWPC1Energy);
    tracks->SetBranchAddress("shortEnergy", shortEnergy);
    tracks->SetBranchAddress("pseudoVertexX", pseudoVertexX);
    tracks->SetBranchAddress("pseudoVertexY", pseudoVertexY);
    tracks->SetBranchAddress("pseudoVertexZ", pseudoVertexZ);
    return true;
}

void TrackInput::GetEntry() {
}
