#ifndef TAGGERINPUT_H
#define TAGGERINPUT_H

#include "InputModule.h"
#include "Rtypes.h"
#include <vector>
#include "TMath.h"

#define GTreeTagger_MAX 4096

namespace ant {
namespace input {

class TaggerInput: public BaseInputModule {
protected:
    Int_t           nTagged = 0;
    Int_t           taggedChannel[GTreeTagger_MAX] = {};
    Double_t        taggedTime[GTreeTagger_MAX] = {};
    Double_t        taggedEnergy[GTreeTagger_MAX] = {};
    Bool_t          hasEnergy = false;
    Double_t        calibration[352] = {};


public:

    TaggerInput();
    virtual ~TaggerInput();

    bool SetupBranches(TreeRequestManager&& input_files);
    void GetEntry();

    Int_t           GetNTagged()                        const	{return nTagged;}
    const	Int_t*          GetTaggedChannel()                  const	{return taggedChannel;}
    Int_t           GetTaggedChannel(const Int_t index) const	{return taggedChannel[index];}
    const	Double_t*       GetTaggedTime()                     const	{return taggedTime;}
    Double_t        GetTaggedTime(const Int_t index)    const	{return taggedTime[index];}
    const	Double_t*       GetTaggedEnergy()                   const	{return taggedEnergy;}
    Double_t        GetTaggedEnergy(const Int_t index)	const	{if(hasEnergy) return taggedEnergy[index]; return calibration[taggedChannel[index]];}
    Bool_t          HasEnergy()                         const   {return hasEnergy;}
    void            SetCalibration(const Int_t nChan, const Double_t *energy);



};
}
}

#endif
