#include "ParticleInput.h"

#include "TTree.h"

using namespace ant;
using namespace input;



ParticleInput::ParticleInput(const std::string name_):
    name(name_)
{

}

ParticleInput::~ParticleInput()
{

}

bool ParticleInput::SetupBranches(TreeRequestManager&& input_files){

    TTree* praticle = input_files.GetTree(name);

    if(praticle==nullptr)
        return false;

    praticle->SetBranchAddress("nParticles",&nParticles);
    praticle->SetBranchAddress("mass", mass);
    praticle->SetBranchAddress("trackIndex", trackIndex);

    return true;
}

void ParticleInput::GetEntry()
{

}
