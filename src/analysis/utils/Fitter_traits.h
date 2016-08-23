#pragma once

#include "APLCON.hpp"
#include "tree/TParticle.h"

namespace ant {
namespace analysis {
namespace utils {

struct Fitter_traits {
    virtual void SetEgammaBeam(double ebeam) =0;
    virtual void SetProton(const TParticlePtr& proton) =0;
    virtual void SetPhotons(const TParticleList& photons) =0;
    virtual void SetZVertexSigma(double sigma) =0;
    virtual bool IsZVertexFitEnabled() const noexcept =0;

    virtual APLCON::Result_t DoFit() =0;

    virtual TParticlePtr GetFittedProton() const =0;
    virtual TParticleList GetFittedPhotons() const =0;
    virtual double GetFittedBeamE() const =0;
    virtual double GetFittedZVertex() const =0;

    virtual ~Fitter_traits() = default;
};

}}} // namespace ant:analysis::utils