#pragma once


#ifndef __CINT__
#include "TParticle.h"
#include <algorithm>
#endif

#include "Rtypes.h"
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
    bool   TouchesHole;

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

    TSimpleParticle() : TLorentzVector() {}
};

}

