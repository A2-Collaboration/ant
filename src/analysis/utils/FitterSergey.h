#pragma once

#include "tree/TParticle.h"

#include "APLCON.hpp"

#include "Rtypes.h"

#include <memory>

namespace ant {
namespace analysis {
namespace utils {

class FitterSergey {

    // use PIMPL idiom to hide all the implementation (copied from Sergey's Acqu)
    class TA2KFitC;
    std::unique_ptr<TA2KFitC> I;

    double Ebeam;
    TParticlePtr Proton;
    TParticleList Photons;

public:
    FitterSergey();
    virtual ~FitterSergey();

    void SetEgammaBeam(double ebeam);
    void SetProton(const TParticlePtr& proton);
    void SetPhotons(const TParticleList& photons);

    APLCON::Result_t DoFit();

    TParticlePtr GetFittedProton() const;
    TParticleList GetFittedPhotons() const;
    double GetFittedBeamE() const;
    double GetFittedZVertex() const;
};



}}} // namespace ant::analysis::utils