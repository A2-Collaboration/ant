#pragma once

#include "InputModule.h"
#include "Rtypes.h"
#include <vector>
#include "TMath.h"


namespace ant {
namespace input {

class DetectorHitInput: public BaseInputModule {
protected:
    Int_t		nNaIHits = 0;
    Int_t		NaIHits[720];
    Int_t		NaICluster[720];
    Int_t		nPIDHits = 0;
    Int_t		PIDHits[24];
    Int_t		nMWPCHits = 0;
    Int_t		MWPCHits[860];
    Int_t		nBaF2Hits = 0;
    Int_t		BaF2Hits[438];
    Int_t		BaF2Cluster[438];
    Int_t		nVetoHits = 0;
    Int_t		VetoHits[438];
public:

    DetectorHitInput();
    virtual ~DetectorHitInput();

    bool SetupBranches(TreeRequestManager&& input_files);
    void GetEntry();

            Int_t		GetNNaIHits()              	const	{return nNaIHits;}
    const	Int_t*		GetNaIHits()           		const	{return NaIHits;}
            Int_t		GetNaIHits(const Int_t index)	const	{return NaIHits[index];}
    const	Int_t*		GetNaICluster()           		const	{return NaICluster;}
            Int_t		GetNaICluster(const Int_t index)	const	{return NaICluster[index];}

            Int_t		GetNPIDHits()      			const	{return nPIDHits;}
    const	Int_t*		GetPIDHits()               	const	{return PIDHits;}
            Int_t		GetPIDHits(const Int_t index)	const	{return PIDHits[index];}

            Int_t		GetNMWPCHits()       			const	{return nMWPCHits;}
    const	Int_t*		GetMWPCHits()                	const	{return MWPCHits;}
            Int_t		GetMWPCHits(const Int_t index)	const	{return MWPCHits[index];}

            Int_t		GetNBaF2Hits()                   const	{return nBaF2Hits;}
    const	Int_t*		GetBaF2Hits()                    const	{return BaF2Hits;}
            Int_t		GetBaF2Hits(const Int_t index)	const	{return BaF2Hits[index];}
    const	Int_t*		GetBaF2Cluster()                    const	{return BaF2Cluster;}
            Int_t		GetBaF2Cluster(const Int_t index)	const	{return BaF2Cluster[index];}

            Int_t		GetNVetoHits()                 const	{return nVetoHits;}
    const	Int_t*		GetVetoHits()                  const	{return VetoHits;}
            Int_t		GetVetoHits(const Int_t index)	const	{return VetoHits[index];}

};
}
}
