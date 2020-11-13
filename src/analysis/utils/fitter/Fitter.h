#pragma once

#include "analysis/utils/Uncertainties.h"
#include "tree/TParticle.h"
#include "base/vec/LorentzVec.h"

#include "APLCON.hpp"

#include <vector>
#include <memory>


namespace ant {
namespace analysis {
namespace utils {

class Fitter {

public:

    static const APLCON::Fit_Settings_t DefaultSettings;

    // Value, Sigma
    struct V_S_t {
        explicit V_S_t(double value = std_ext::NaN) : Value(value) {}
        double Value;
        double Sigma = std_ext::NaN;
        void SetValueSigma(double value, double sigma) {
            Value = value;
            Sigma = sigma;
        }
        operator const double&() const noexcept {
            return Value;
        }
    };

    // Value, Sigma, Pull
    struct V_S_P_t : V_S_t {
        using V_S_t::V_S_t; // use base class ctor
        double Pull = std_ext::NaN;

        template<std::size_t N>
        std::tuple<double&> linkFitter() noexcept {
            // the following get<N> assumes this order of indices
            static_assert(APLCON::ValueIdx==0,"");
            static_assert(APLCON::SigmaIdx==1,"");
            static_assert(APLCON::PullIdx ==2,"");
            // the extra std::tie around std::get is for older compilers...
            return std::tie(std::get<N>(std::tie(Value, Sigma, Pull)));
        }
    };

    struct Z_Vertex_t : V_S_P_t {
        explicit Z_Vertex_t(bool isEnabled) :
            V_S_P_t(0), // always init Z_Vertex to 0, as it's used by FitParticle::GetLorentzVec
            IsEnabled(isEnabled) {}
        const bool IsEnabled;  /// \todo user might want to change it from fit to fit
        double Sigma_before = std_ext::NaN;

        template<size_t innerIdx>
        APLCON::Variable_Settings_t getFitterSettings(size_t outerIdx) const noexcept {
            static_assert(innerIdx==0,""); // just be sure
            (void)outerIdx; // unused, provided to user method for completeness
            APLCON::Variable_Settings_t settings;
            // not enabled z-vertex is fixed at zero
            if(!IsEnabled)
                settings.StepSize = 0;
            return settings;
        }
    private:
        using V_S_P_t::operator const double&;
    };

    struct FitParticle
    {
        TParticlePtr Particle; // pointer to unfitted particle

        TParticlePtr AsFitted() const;

        double GetShowerDepth() const { return ShowerDepth; }
        std::vector<double> GetValues_before() const;
        std::vector<double> GetSigmas_before() const;
        std::vector<double> GetPulls() const;

        template<std::size_t N>
        std::tuple<double&, double&, double&, double&> linkFitter() noexcept {
            if(N==APLCON::ValueIdx)
                return std::tie(Vars[0].Value, Vars[1].Value, Vars[2].Value, Vars[3].Value);
            else if(N==APLCON::SigmaIdx)
                return std::tie(Vars[0].Sigma, Vars[1].Sigma, Vars[2].Sigma, Vars[3].Sigma);
            else // N == APLCON::PullIdx
                return std::tie(Pulls[0], Pulls[1], Pulls[2], Pulls[3]);
        }

    protected:

        friend class KinFitter;
        friend class TreeFitter;
        friend class SigmaFitter;
        friend class NoProtonFitter;

        void Set(const TParticlePtr& p, const UncertaintyModel& model);
        void SetFittedZVertex(double zvertex) { Fitted_Z_Vertex = zvertex; }
        ant::LorentzVec GetLorentzVec(double zvertex) const noexcept;

        void SetEk(double Ek) noexcept { Vars[0].Value = 1.0/Ek; }
        bool IsEkUnmeasured() const noexcept { return Vars[0].Sigma == 0; }

        double ShowerDepth = std_ext::NaN;
        double Fitted_Z_Vertex = std_ext::NaN;

        // CB and TAPS both use 4 parameters
        // (but with different meaning, see Set() method and GetLorentzVec())
        std::array<V_S_t,  4> Vars;
        std::array<V_S_t,  4> Vars_before;
        std::array<double, 4> Pulls;
    };

    struct Exception : std::runtime_error {
        using std::runtime_error::runtime_error;
    };

protected:

    // one should derive from that class
    Fitter() = default;


private:
    static APLCON::Fit_Settings_t MakeDefaultSettings();

};

}}} // namespace ant::analysis::utils


