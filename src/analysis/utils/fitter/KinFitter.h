#pragma once

#include "Fitter.h"

namespace ant {
namespace analysis {
namespace utils {

class KinFitter : public Fitter
{
public:

    /**
     * @brief KinFitter applies energy-momentum constraint to proton/photons using incoming beam
     * @param uncertainty_model model to obtain uncertainties
     * @param settings tune the underlying APLCON fitter
     */
    explicit KinFitter(UncertaintyModelPtr uncertainty_model = nullptr,
                       bool fit_Z_vertex = false,
                       const APLCON::Fit_Settings_t& settings = DefaultSettings
             );

    void SetZVertexSigma(double sigma);
    bool IsZVertexFitEnabled() const noexcept;
    bool IsZVertexUnmeasured() const;

    TParticlePtr  GetFittedProton() const;
    TParticleList GetFittedPhotons() const;
    double GetFittedBeamE() const;
    TParticlePtr GetFittedBeamParticle() const;
    double GetFittedZVertex() const;

    double GetBeamEPull() const;
    double GetZVertexPull() const;

    std::vector<FitParticle> GetFitParticles() const;

    APLCON::Result_t DoFit(double ebeam, const TParticlePtr& proton, const TParticleList& photons);

    void SetUncertaintyModel(const UncertaintyModelPtr& uncertainty_model) {
        Model = uncertainty_model;
    }

protected:

    void PrepareFit(double ebeam,
                    const TParticlePtr& proton,
                    const TParticleList& photons);


    struct BeamE_t : V_S_P_t {
        double Value_before = std_ext::NaN;
        LorentzVec GetLorentzVec() const noexcept;
        void SetValueSigma(double value, double sigma) {
            V_S_P_t::SetValueSigma(value, sigma);
            Value_before = Value;
        }
    };

    using Proton_t = FitParticle;
    using Photons_t = std::vector<FitParticle>;

    BeamE_t    BeamE;
    Proton_t   Proton;
    Photons_t  Photons;
    Z_Vertex_t Z_Vertex;

    APLCON::Fitter<BeamE_t, Proton_t, Photons_t, Z_Vertex_t> aplcon;

    // make constraint a static function, then we can use the typedefs
    static std::array<double, 4> constraintEnergyMomentum(const BeamE_t& beam, const Proton_t& proton,
                                                          const Photons_t& photons, const Z_Vertex_t&);


private:
    UncertaintyModelPtr Model;
};

}}} // namespace ant::analysis::utils