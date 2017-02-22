#pragma once

#include "analysis/utils/Uncertainties.h"
#include "tree/TParticle.h"
#include "base/vec/LorentzVec.h"

#include "APLCON.hpp"

namespace ant {
namespace analysis {
namespace utils {

class Fitter {

public:

    struct Exception : std::runtime_error {
        using std::runtime_error::runtime_error;
    };

    static const APLCON::Fit_Settings_t DefaultSettings;

    Fitter(const Fitter&) = delete;
    Fitter& operator=(const Fitter&) = delete;
    virtual ~Fitter();

    struct FitVariable {
        // initialize linked values at zero, important for Z Vertex
        double Value = 0;
        double Value_before = std_ext::NaN;
        double Sigma = 0;
        double Sigma_before = std_ext::NaN;
        double Pull = 0;
        void SetValueSigma(double v, double s) {
            Value = v;
            Value_before = v;
            Sigma = s;
            Sigma_before = s;
        }
    };

    struct FitParticle
    {
        TParticlePtr Particle; // pointer to unfitted particle

        std::vector<FitVariable> Vars;

        TParticlePtr AsFitted() const;

        double GetShowerDepth() const;
        std::vector<double> GetValues() const;
        std::vector<double> GetSigmas() const;
        std::vector<double> GetPulls() const;

        FitParticle(const std::string& name,
                    APLCON& aplcon,
                    std::shared_ptr<FitVariable> z_vertex);
        virtual ~FitParticle();

    protected:

        friend class Fitter;
        friend class KinFitter;
        friend class TreeFitter;

        void Set(const TParticlePtr& p, const UncertaintyModel& uncertainty);

        const std::string Name;
        const std::shared_ptr<const FitVariable> Z_Vertex;

        ant::LorentzVec GetLorentzVec(const std::vector<double>& values,
                                      double z_vertex) const;

    private:
        double ShowerDepth = std_ext::NaN;
    };

protected:

    Fitter(const std::string& fittername,
           const APLCON::Fit_Settings_t& settings,
           UncertaintyModelPtr uncertainty_model);

    Fitter(Fitter&&) = default;
    Fitter& operator=(Fitter&&) = default;

    UncertaintyModelPtr uncertainty;
    std::unique_ptr<APLCON> aplcon;

private:
    static APLCON::Fit_Settings_t MakeDefaultSettings();

};

}}} // namespace ant::analysis::utils


