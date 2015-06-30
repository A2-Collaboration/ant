#include "PlutoInput.h"
#include "TTree.h"
#include "PParticle.h"
#include "TClonesArray.h"

using namespace ant;
using namespace ant::input;

bool PlutoInput::SetupBranches(TreeRequestManager&& input_files) {

    TTree* data = input_files.GetTree("data");
    if(data==nullptr)
        return false;

    data->SetBranchAddress("Particles",         &PlutoMCTrue);
    data->SetBranchAddress("plutoID",           &plutoID);
    data->SetBranchAddress("plutoRandomID", 	&plutoRandomID);

    return true;
}

void PlutoInput::GetEntry()
{
    const Long64_t entries = PlutoMCTrue->GetEntries();
    particles.clear();
    particles.reserve(size_t(entries));

    for( Long64_t i=0; i < entries; ++i) {

        const PParticle* const particle = dynamic_cast<const PParticle*>((*PlutoMCTrue)[i]);

        if(particle) {
            particles.push_back(particle);
        }
    }

}


PlutoInput::PlutoInput()
{

}

PlutoInput::~PlutoInput()
{

}
