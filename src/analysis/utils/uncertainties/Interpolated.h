#pragma once

#include "analysis/utils/Uncertainties.h"
#include "base/ClippedInterpolatorWrapper.h"

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
struct Interpolated : public UncertaintyModel {
public:

    Interpolated(UncertaintyModelPtr starting_uncertainty_, bool use_proton_sigmaE_ = false);
    virtual ~Interpolated();

    Uncertainties_t GetSigmas(const TParticle &particle) const override;

    void LoadSigmas(const std::string& filename);

    bool HasLoadedSigmas() const {
        return loaded_sigmas;
    }

    enum class Type_t {
        Data, MC
    };

    static std::shared_ptr<Interpolated> makeAndLoad(
            UncertaintyModelPtr default_model = nullptr,
            bool use_proton_sigmaE = false) {
        return makeAndLoad(Type_t::Data, default_model, use_proton_sigmaE);
    }

    static std::shared_ptr<Interpolated> makeAndLoad(
            Type_t type,
            UncertaintyModelPtr default_model = nullptr,
            bool use_proton_sigmaE = false);

    friend std::ostream& operator<<(std::ostream& stream, const Interpolated& o);

protected:
    const UncertaintyModelPtr starting_uncertainty;
    const bool use_proton_sigmaE;

    bool loaded_sigmas = false;

    static std::unique_ptr<const Interpolator2D> LoadInterpolator(const WrapTFile& file, const std::string& prefix);

    struct EkThetaPhiR {

        ClippedInterpolatorWrapper Ek;
        ClippedInterpolatorWrapper Theta;
        ClippedInterpolatorWrapper Phi;
        ClippedInterpolatorWrapper CB_R;
        ClippedInterpolatorWrapper ShowerDepth;

        void SetUncertainties(Uncertainties_t& u, const TParticle& particle) const;
        void Load(const WrapTFile& file, const std::string& prefix);
    };
    friend std::ostream& operator<<(std::ostream& stream, const EkThetaPhiR& o);

    struct EkRxyPhiL {

        ClippedInterpolatorWrapper Ek;
        ClippedInterpolatorWrapper TAPS_Rxy;
        ClippedInterpolatorWrapper Phi;
        ClippedInterpolatorWrapper TAPS_L;
        ClippedInterpolatorWrapper ShowerDepth;

        void SetUncertainties(Uncertainties_t& u, const TParticle& particle) const;
        void Load(const WrapTFile& file, const std::string& prefix);
    };
    friend std::ostream& operator<<(std::ostream& stream, const EkRxyPhiL& o);

    EkThetaPhiR cb_photon;
    EkRxyPhiL   taps_photon;
    EkThetaPhiR cb_proton;
    EkRxyPhiL   taps_proton;

};

}}}} // namespace ant::analysis::utils::UncertaintyModels
