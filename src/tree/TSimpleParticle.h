#pragma once


#include "base/types.h"

#include "base/Tree.h"

#include "TParticle.h"
#include "base/vec/LorentzVec.h"




namespace ant {

/**
 * @brief Base TParticle class
 */
struct TSimpleParticle : LorentzVec {

    double VetoE;
    double ShortE;

    int    ClusterSize;

    TSimpleParticle(const TParticle& particle);

    TSimpleParticle(const LorentzVec& lorentzvector,
              const double vetoE, const double shortE,
              const int clusterSize):
        LorentzVec(lorentzvector),
        VetoE(vetoE), ShortE(shortE), ClusterSize(clusterSize){}

    friend std::ostream& operator<<(std::ostream& stream, const TSimpleParticle& o);

    TSimpleParticle(const TSimpleParticle&) = default;
    TSimpleParticle& operator= (const TSimpleParticle&) = default;
    TSimpleParticle(TSimpleParticle&&) = default;
    TSimpleParticle& operator= (TSimpleParticle&&) = default;

    TSimpleParticle() : LorentzVec(),
        VetoE(std_ext::NaN), ShortE(std_ext::NaN), ClusterSize(){}

};

}

