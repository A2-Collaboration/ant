#pragma once

#include "tree/TParticle.h"
#include "base/Detector_t.h"

namespace ant {
namespace analysis {
namespace utils {

/**
 * @brief Uncertainties for E, theta, and phi
 * and some CB/TAPS specific values
 */
struct Uncertainties_t {

    double sigmaEk;
    double sigmaTheta;
    double sigmaPhi;

    Detector_t::Any_t Detector;
    double ShowerDepth;
    double sigmaCB_R;
    double sigmaTAPS_Rxy;
    double sigmaTAPS_L;

    Uncertainties_t(double Ek = std_ext::NaN,
                    double Theta = std_ext::NaN,
                    double Phi = std_ext::NaN,
                    Detector_t::Any_t detector = Detector_t::Any_t::None,
                    double showerDepth = std_ext::NaN,
                    double CB_R = std_ext::NaN,
                    double TAPS_Rxy = std_ext::NaN,
                    double TAPS_L = std_ext::NaN) :
        sigmaEk(Ek), sigmaTheta(Theta), sigmaPhi(Phi),
        Detector(detector), ShowerDepth(showerDepth),
        sigmaCB_R(CB_R), sigmaTAPS_Rxy(TAPS_Rxy), sigmaTAPS_L(TAPS_L)
    {}
};

/**
 * @brief Virtual base class for different Uncertainty Models for fitter.
 *        Implement at least the GetSigmas() method.
 * @see UncertaintyModels::Constant
 * @see UncertaintyModels::MCExtracted
 * @see UncertaintyModels::Interpolated
 */
class UncertaintyModel {
public:
    UncertaintyModel();
    virtual ~UncertaintyModel();
    virtual Uncertainties_t GetSigmas(const TParticle& particle) const =0;
    virtual double GetBeamEnergySigma(double photon_energy) const;
    struct Exception : std::runtime_error {
        using std::runtime_error::runtime_error;
    };
protected:
    std::shared_ptr<TaggerDetector_t> tagger;
};
using UncertaintyModelPtr = std::shared_ptr<const UncertaintyModel>;

}}} // namespace ant::analysis::utils
