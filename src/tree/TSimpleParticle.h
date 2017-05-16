#pragma once


#ifndef __CINT__
#include "TParticle.h"
#include <algorithm>
#endif

#include "TLorentzVector.h"



namespace ant {

/**
 * @brief Base TSimpleParticle class
 */
struct TSimpleParticle : TLorentzVector {
    ClassDef(TSimpleParticle,1)

    double Time;

    double VetoE;
    double TrackerE;

    double ShortE;

    double Mass;

    int    ClusterSize;

    double Ek() const { return E() - Mass; }

#ifndef __CINT__
    TSimpleParticle(const TParticle& particle);
    static std::vector<TSimpleParticle> TransformParticleList(const TParticleList& particles)
    {
        std::vector<TSimpleParticle> vpart(particles.size());
        std::transform(particles.begin(),particles.end(),vpart.begin(),
                       [](const TParticlePtr& ph){return TSimpleParticle(*ph);});
        return vpart;
    }
#endif

    TSimpleParticle(const TLorentzVector& lorentzvector,
                    const double time,
                    const double vetoE, const double trackerE,
                    const double shortE,
                    const double mass,
                    const int clusterSize):
        TLorentzVector(lorentzvector),
        Time(time),
        VetoE(vetoE), TrackerE(trackerE),
        ShortE(shortE),
        Mass(mass),
        ClusterSize(clusterSize){}

    TSimpleParticle(): TLorentzVector(){}


};

}

