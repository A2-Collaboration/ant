#pragma once


#ifndef __CINT__
#include "TParticle.h"
#endif

#include "TLorentzVector.h"



namespace ant {

/**
 * @brief Base TParticle class
 */
struct TSimpleParticle : TLorentzVector {
    ClassDef(TSimpleParticle,1)

    double Time;

    double VetoE;
    double ShortE;

    double Mass;

    int    ClusterSize;

#ifndef __CINT__
    TSimpleParticle(const TParticle& particle);
#endif

    TSimpleParticle(const TLorentzVector& lorentzvector,
                    const double time,
                    const double vetoE, const double shortE,
                    const double mass,
                    const int clusterSize):
        TLorentzVector(lorentzvector),
        Time(time),
        VetoE(vetoE), ShortE(shortE),
        Mass(mass),
        ClusterSize(clusterSize){}

    TSimpleParticle(): TLorentzVector(){}


};

}

