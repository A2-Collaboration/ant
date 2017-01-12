#pragma once

#include "analysis/utils/Uncertainties.h"

namespace ant {

class Interpolator2D;
class WrapTFile;

namespace analysis {
namespace utils {
namespace UncertaintyModels {

/**
 * @brief Uncertainties with interpolated surfaces in (E,theta) plane,
 * determined with iterative procedure
 * @see progs/Ant-makeSigmas.cc
 * @see src/analysis/physics/common/InterpolatedPulls.h
 */
struct Interpolated : public UncertaintyModel, public ant::printable_traits {
public:

    Interpolated(UncertaintyModelPtr starting_uncertainty_, bool use_proton_sigmaE_ = false);
    virtual ~Interpolated();

    Uncertainties_t GetSigmas(const TParticle &particle) const override;

    void LoadSigmas(const std::string& filename);

    bool HasLoadedSigmas() const {
        return loaded_sigmas;
    }

    static std::shared_ptr<Interpolated> makeAndLoad(UncertaintyModelPtr default_model = nullptr,
                                                     bool use_proton_sigmaE = false);

    struct ClippedInterpolatorWrapper : ant::printable_traits {
        using interpolator_ptr_t = std::unique_ptr<const Interpolator2D>;

        interpolator_ptr_t interp;

        struct boundsCheck_t : ant::printable_traits {
            interval<double> range;
            mutable unsigned underflow = 0;
            mutable unsigned unclipped = 0;
            mutable unsigned overflow  = 0;

            double clip(double v) const;

            boundsCheck_t(const interval<double> r): range(r) {}

            std::ostream& Print(std::ostream& stream) const override;
        };

        boundsCheck_t xrange;
        boundsCheck_t yrange;

        ClippedInterpolatorWrapper(interpolator_ptr_t i);
        ClippedInterpolatorWrapper();
        ~ClippedInterpolatorWrapper();
        double GetPoint(double x, double y) const;

        void setInterpolator(interpolator_ptr_t i);

        std::ostream& Print(std::ostream& stream) const override;

    };

    std::ostream& Print(std::ostream& stream) const override;

protected:
    const UncertaintyModelPtr starting_uncertainty;
    const bool use_proton_sigmaE;

    bool loaded_sigmas = false;

    static std::unique_ptr<const Interpolator2D> LoadInterpolator(const WrapTFile& file, const std::string& prefix);

    struct EkThetaPhiR : ant::printable_traits {

        ClippedInterpolatorWrapper Ek;
        ClippedInterpolatorWrapper Theta;
        ClippedInterpolatorWrapper Phi;
        ClippedInterpolatorWrapper CB_R;
        ClippedInterpolatorWrapper ShowerDepth;

        void SetUncertainties(Uncertainties_t& u, const TParticle& particle) const;
        void Load(const WrapTFile& file, const std::string& prefix);
        std::ostream& Print(std::ostream& stream) const override;
    };

    struct EkRxyPhiL : ant::printable_traits {

        ClippedInterpolatorWrapper Ek;
        ClippedInterpolatorWrapper TAPS_Rxy;
        ClippedInterpolatorWrapper Phi;
        ClippedInterpolatorWrapper TAPS_L;
        ClippedInterpolatorWrapper ShowerDepth;

        void SetUncertainties(Uncertainties_t& u, const TParticle& particle) const;
        void Load(const WrapTFile& file, const std::string& prefix);
        std::ostream& Print(std::ostream& stream) const override;
    };

    EkThetaPhiR cb_photon;
    EkRxyPhiL   taps_photon;
    EkThetaPhiR cb_proton;
    EkRxyPhiL   taps_proton;
};

}}}} // namespace ant::analysis::utils::UncertaintyModels