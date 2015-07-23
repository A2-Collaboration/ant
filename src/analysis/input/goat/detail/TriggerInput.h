#pragma once

#include "InputModule.h"
#include "Rtypes.h"
#include <vector>
#include "TMath.h"

#define GTreeTrigger_MAX 256

namespace ant {
namespace input {

class TriggerInput: public BaseInputModule {
protected:
    Double_t 	energySum   = 0.0;
    Int_t 		multiplicity = 0;
    Int_t 		nTriggerPattern = 0;
    Int_t 		triggerPattern[GTreeTrigger_MAX];
    Int_t 		nErrors = 0;
    Int_t 		errorModuleID[GTreeTrigger_MAX];
    Int_t 		errorModuleIndex[GTreeTrigger_MAX];
    Int_t 		errorCode[GTreeTrigger_MAX];
    Bool_t      helicity = 0;
    Long64_t    MC_evt_id = -1;
    Long64_t    MC_rnd_id = -1;
    Bool_t      hasHelicity = false;
    Bool_t      hasMCID = false;

public:

    TriggerInput();
    virtual ~TriggerInput();

    bool SetupBranches(TreeRequestManager&& input_files);
    void GetEntry();

    virtual void        Clear()                       {nTriggerPattern = 0; nErrors = 0;}
            Int_t 		GetMultiplicity()       const {return multiplicity;}
            Double_t	GetEnergySum()          const {return energySum;}
            Int_t		GetNTriggerPattern()    const {return nTriggerPattern;}
    const   Int_t*		GetTriggerPattern()     const {return triggerPattern;}
            Int_t		GetNErrors()            const {return nErrors;}
    const   Int_t*		GetErrorModuleID()      const {return errorModuleID;}
    const   Int_t*		GetErrorModuleIndex()   const {return errorModuleIndex;}
    const   Int_t*		GetErrorCode()          const {return errorCode;}
            Bool_t 	    GetHelicity()    	    const {return helicity;}
            Bool_t 	    HasHelicity()    	    const {return hasHelicity;}
            Long64_t    GetMCTrueEventID()      const {return MC_evt_id;}
            Long64_t    GetMCRandomID()         const {return MC_rnd_id;}
};
}
}
