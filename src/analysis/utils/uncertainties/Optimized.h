#pragma once

#include "analysis/utils/Uncertainties.h"

namespace ant {
namespace analysis {
namespace utils {
namespace UncertaintyModels {

/**
 * @brief Optimized / adjusted to data uncertainty model:
 *
 * Uses fixed/constant values for protons: Energy always 0 -> unmeasured,
 * different values for theta and phi in CB and TAPS.
 *
 * Photons:
 * CB:
 *  simga E = A*E / (E[GeV])**b
 *
 *  sigma Theta = constant + A*sin(theta)
 *
 *  sigma Phi   = sPhi / sin(theta)
 *
 * TAPS:
 *   constant values for now.
 *
 */
struct Optimized : UncertaintyModel {

    double cb_photon_theta_const = 0.0;
    double cb_photon_theta_Sin   = 0.0;
    double cb_photon_phi         = 0.0;

    double cb_photon_E_rel = 0.0;
    double cb_photon_E_exp = 0.0;
    double cb_photon_E_lin = 0.0;

    Uncertainties_t cb_proton = {};


    double taps_photon_E_rel = 0.0;
    double taps_photon_E_exp = 0.0;
    double taps_photon_E_lin = 0.0;
    double taps_photon_theta = 0.0;
    double taps_photon_phi   = 0.0;

    Uncertainties_t taps_proton = {};

    Optimized();

    virtual ~Optimized();

    Uncertainties_t GetSigmas(const TParticle& particle) const;
    static double dThetaSin(const double theta, const double offset, const double thetapart) noexcept;

    static double dE(const double E, const double rel, const double exp, const double reloffset) noexcept;

    std::string to_string_simple() const;

    /**
     * @brief serialize to JSON string
     * @return
     */
    std::string to_string() const;

    /**
     * @brief serialzie to JSON string without whitespace
     * @return
     */
    std::string to_string_short() const;

    /**
     * @brief load from a JSON string
     * @param data
     */
    void load_from_string(const std::string& data);

    void load_from_string_simple(const std::string& data);

    bool operator==(const Optimized& other) const noexcept;
    bool operator!=(const Optimized& other) const noexcept;

protected:
    const static std::string separator;

    void ReadToken(const std::string& token);
};

/**
 * @brief Values for the optimized model, obtainded from KinFitPi0 on data
 */
struct Optimized_Oli1 : Optimized {
    Optimized_Oli1(double relative_scale = 1.0, bool use_measured_proton_TAPS = false);
    virtual ~Optimized_Oli1();
};

}}}} // namespace ant::analysis::utils::UncertaintyModels