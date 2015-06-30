#ifndef TRACKINPUT_H
#define TRACKINPUT_H

#include "FileManager.h"
#include "InputModule.h"
#include "Rtypes.h"
#include <vector>
#include "TMath.h"

#define GTreeTrack_MAX 128

namespace ant {
namespace input {

class TrackInput: public BaseInputModule {
public:
    enum
    {
        DETECTOR_NONE = 0,
        DETECTOR_NaI = 1,
        DETECTOR_PID = 2,
        DETECTOR_MWPC = 4,
        DETECTOR_BaF2 = 8,
        DETECTOR_PbWO4 = 16,
        DETECTOR_Veto = 32,
    };
protected:
    Int_t		nTracks = 0;
    Double_t	clusterEnergy[GTreeTrack_MAX];
    Double_t 	theta[GTreeTrack_MAX];
    Double_t 	phi[GTreeTrack_MAX];
    Double_t	time[GTreeTrack_MAX];
    Int_t       clusterSize[GTreeTrack_MAX];
    Int_t       centralCrystal[GTreeTrack_MAX];
    Int_t       centralVeto[GTreeTrack_MAX];
    Int_t       detectors[GTreeTrack_MAX];
    //Charged detector energies
    Double_t	vetoEnergy[GTreeTrack_MAX];
    Double_t	MWPC0Energy[GTreeTrack_MAX];
    Double_t	MWPC1Energy[GTreeTrack_MAX];
    //TAPS PSA Short-gate Energy
    Double_t	shortEnergy[GTreeTrack_MAX];
    //Pseudo vertex information
    Double_t    pseudoVertexX[GTreeTrack_MAX];
    Double_t    pseudoVertexY[GTreeTrack_MAX];
    Double_t    pseudoVertexZ[GTreeTrack_MAX];

public:

    TrackInput();
    virtual ~TrackInput();

    bool SetupBranches(TreeRequestManager&& input_files);
    void GetEntry();

            Int_t           GetDetectors(const Int_t index)     const	{return detectors[index];}
            Int_t           GetClusterSize(const Int_t index)   const 	{return clusterSize[index];}
            Int_t           GetCentralCrystal(const Int_t index)  const {return centralCrystal[index];}
            Int_t           GetCentralVeto(const Int_t index)   const 	{return centralVeto[index];}
            Double_t        GetVetoEnergy(const Int_t index)    const	{return vetoEnergy[index];}
            Double_t        GetClusterEnergy(const Int_t index) const	{return clusterEnergy[index];}
            Int_t           GetNTracks()                        const	{return nTracks;}
    inline  Int_t           GetNCB()                            const;
    inline  Int_t           GetNTAPS()                          const;
    inline  Bool_t          HasCB(const Int_t index)            const;
    inline  Bool_t          HasTAPS(const Int_t index)          const;
            Double_t        GetPhi(const Int_t index)        const	{return phi[index] * TMath::DegToRad();}
            Double_t        GetTheta(const Int_t index)      const	{return theta[index] * TMath::DegToRad();}
            Double_t        GetTime(const Int_t index)          const	{return time[index];}
            Double_t        GetMWPC0Energy(const Int_t index)         const	{return MWPC0Energy[index];}
            Double_t        GetMWPC1Energy(const Int_t index)         const	{return MWPC1Energy[index];}
            Double_t        GetShortEnergy(const Int_t index)         const	{return shortEnergy[index];}
            Double_t        GetPseudoVertexX(const Int_t index)       const	{return pseudoVertexX[index];}
            Double_t        GetPseudoVertexY(const Int_t index)       const	{return pseudoVertexY[index];}
            Double_t        GetPseudoVertexZ(const Int_t index)       const	{return pseudoVertexZ[index];}
};
}
}

#endif
