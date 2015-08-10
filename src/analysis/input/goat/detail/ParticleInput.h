#pragma once

#include "InputModule.h"
#include "Rtypes.h"
#include <vector>
#include "TMath.h"

#define GTreeTrack_MAX 128

namespace ant {
namespace analysis {
namespace input {

class ParticleInput: public BaseInputModule {
protected:
    std::string         name;
    Int_t               nParticles = 0;
    Double_t            mass[GTreeTrack_MAX] = {};
    Int_t               trackIndex[GTreeTrack_MAX] = {};  // index of the corresponding tack in the track list, -1 => No track

public:

    ParticleInput(const std::string name_);
    virtual ~ParticleInput();

    bool SetupBranches(TreeRequestManager&& input_files);
    void GetEntry();

    Int_t           GetNParticles()                     const	{return nParticles;}
    Double_t        GetMass(const Int_t index)          const   {return mass[index];}
    Int_t           GeTCandidateIndex(const Int_t index)    const   {return trackIndex[index];}

};
}
}
}
