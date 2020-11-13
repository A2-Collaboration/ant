#pragma once

#include "Fitter.h"

namespace ant {
namespace analysis {
namespace utils {

class NoProtonFitter : public Fitter
{
public:

    /**
     * @brief NoProtonFitter applies energy-momentum constraint to proton/photons where proton is set as unmeasured, STILL NOT COMPLETELY TESTED
     * @param uncertainty_model model to obtain uncertainties
     * @param settings tune the underlying APLCON fitter
     */
    explicit NoProtonFitter(UncertaintyModelPtr uncertainty_model = nullptr,
                       bool fit_Z_vertex = false,
                       const APLCON::Fit_Settings_t& settings = DefaultSettings
             );

    void SetZVertexSigma(double sigma);
    bool IsZVertexFitEnabled() const noexcept;
    bool IsZVertexUnmeasured() const;
    void SetTarget(double length, double center = 0.);

    LorentzVec  GetFittedProton() const;
    TParticleList GetFittedPhotons() const;
    double GetFittedBeamE() const;
    TParticlePtr GetFittedBeamParticle() const;
    double GetFittedZVertex() const;

    double GetBeamEPull() const;
    double GetZVertexPull() const;

    std::vector<FitParticle> GetFitParticles() const;

    APLCON::Result_t DoFit(double ebeam, const TParticleList& photons);

    void SetUncertaintyModel(const UncertaintyModelPtr& uncertainty_model) {
        Model = uncertainty_model;
    }

protected:

    void PrepareFit(double ebeam,
                    const TParticleList& photons);

    struct BeamE_t : V_S_P_t {
        double Value_before = std_ext::NaN;
        LorentzVec GetLorentzVec() const noexcept;
        void SetValueSigma(double value, double sigma) {
            V_S_P_t::SetValueSigma(value, sigma);
            Value_before = Value;
        }
    };

    struct Target_t {
        double length = std_ext::NaN;
        double center = std_ext::NaN;

        double start() const { return center - length/2.; }
        double end() const { return center + length/2.; }
    };

    struct FitProtonUnMeas
    {
        ant::LorentzVec AsFitted() const;

        std::vector<double> GetValues_before() const;
        std::vector<double> GetSigmas_before() const;
        std::vector<double> GetPulls() const;

        template<std::size_t N>
        std::tuple<double&, double&, double&> linkFitter() noexcept {
            if(N==APLCON::ValueIdx)
                return std::tie(Vars[0].Value, Vars[1].Value, Vars[2].Value);
            else if(N==APLCON::SigmaIdx)
                return std::tie(Vars[0].Sigma, Vars[1].Sigma, Vars[2].Sigma);
            else // N == APLCON::PullIdx
                return std::tie(Pulls[0], Pulls[1], Pulls[2]);
        }

    protected:

        friend class NoProtonFitter;

        void Set(const LorentzVec);
        ant::LorentzVec GetLorentzVec() const noexcept;

        // For an unmeasured particle, 3 "true"-parameters suffices
        std::array<V_S_t,  3> Vars;
        std::array<V_S_t,  3> Vars_before;
        std::array<double, 3> Pulls;
    };

    struct Exception : std::runtime_error {
                using std::runtime_error::runtime_error;
            };

    using Photons_t = std::vector<FitParticle>;
    using Proton_t = FitProtonUnMeas;

    BeamE_t BeamE;
    Proton_t Proton;
    Photons_t Photons;
    Z_Vertex_t Z_Vertex;
    Target_t Target;

    APLCON::Fitter<BeamE_t, Proton_t, Photons_t, Z_Vertex_t> aplcon;

    // make constraint a static function, then we can use the typedefs
    static std::array<double, 4> constraintEnergyMomentum(const BeamE_t& beam, const Proton_t& proton,
                                                          const Photons_t& photons, const Z_Vertex_t&);

private:
    UncertaintyModelPtr Model;
};

}}} // namespace ant::analysis::utils
