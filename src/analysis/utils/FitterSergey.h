#pragma once

#include "Fitter_traits.h"

#include "tree/TParticle.h"

#include "APLCON.hpp"

#include "Rtypes.h"

#include <memory>

namespace ant {
namespace analysis {
namespace utils {

class FitterSergey : public Fitter_traits {

    // use PIMPL idiom to hide all the implementation (copied from Sergey's Acqu)
    class TA2KFitC;
    std::unique_ptr<TA2KFitC> I;

    double Ebeam;
    TParticlePtr Proton;
    TParticleList Photons;

public:
    FitterSergey();
    virtual ~FitterSergey();

    virtual void SetEgammaBeam(double ebeam) override;
    virtual void SetProton(const TParticlePtr& proton) override;
    virtual void SetPhotons(const TParticleList& photons) override;

    virtual APLCON::Result_t DoFit() override;

    virtual TParticlePtr GetFittedProton() const override;
    virtual TParticleList GetFittedPhotons() const override;
    virtual double GetFittedBeamE() const override;
    virtual double GetFittedZVertex() const override;
};



}}} // namespace ant::analysis::utils