#pragma once

#include "Fitter.h"

namespace ant {
namespace analysis {
namespace utils {

class KinFitter : public Fitter
{
public:

    /**
     * @brief KinFitter::KinFitter
     * @param name Name for the fitter
     * @param numGammas number of photons involved in the fit
     * @param Uncertainty_model model to predict uncertainties
     * @param settings
     */
    KinFitter(const std::string& name,
              unsigned numGammas,
              UncertaintyModelPtr uncertainty_model,
              bool fit_Z_vertex = false,
              const APLCON::Fit_Settings_t& settings = DefaultSettings
              );

    virtual ~KinFitter();

    KinFitter(const KinFitter&) = delete;
    KinFitter& operator=(const KinFitter&) = delete;
    KinFitter(KinFitter&&) = default;
    KinFitter& operator=(KinFitter&&) = default;

    void SetEgammaBeam(double ebeam);
    void SetProton(const TParticlePtr& proton);
    void SetPhotons(const TParticleList& photons);

    void SetZVertexSigma(double sigma);
    bool IsZVertexFitEnabled() const noexcept;

    TParticlePtr GetFittedProton() const;
    TParticleList GetFittedPhotons() const;
    double GetFittedBeamE() const;
    TParticlePtr GetFittedBeamParticle() const;
    double GetFittedZVertex() const;

    double GetBeamEPull() const;
    double GetZVertexPull() const;

    std::vector<double> GetProtonPulls() const;
    /**
     * @brief GetPhotonsPulls
     * @return matrix with first index specifying parameter (0...3), second the photons.
     * in congruence with GetProtonPulls
     */
    std::vector<std::vector<double>> GetPhotonsPulls() const;

    /**
     * @brief GetFitParticles returns as first item the proton, then all n photons
     * @return FitParticle contain all info about the fitted state
     */
    std::vector<FitParticle> GetFitParticles() const;

    APLCON::Result_t DoFit();

protected:

    struct BeamE_t : FitVariable {
        const std::string Name = "Beam";
    };

    struct Z_Vertex_t : FitVariable {
        const std::string Name = "ZVertex";
    };

    // it's pretty important that those things are pointers,
    // since the members are linked to APLCON in ctor!
    // A move/copy of those members may not happen, so we just
    // point to their fixed location in memory.
    std::vector<std::shared_ptr<FitParticle>> Photons;
    std::shared_ptr<FitParticle> Proton;
    std::unique_ptr<BeamE_t>    BeamE;
    std::shared_ptr<Z_Vertex_t> Z_Vertex;

    static LorentzVec MakeBeamLorentzVec(double BeamE);

};

}}} // namespace ant::analysis::utils