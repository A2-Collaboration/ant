#pragma once

#include "Fitter.h"

namespace ant {

class WrapTFile;

namespace analysis {
namespace utils {

namespace UncertaintyModels {

/**
 * @brief Constant/static kin fitter uncertainty model. has a fixed value for {cb,taps}{photon,proton}{E,theta,phi}
 */
struct Constant : public Fitter::UncertaintyModel {
public:

    struct Exception : std::runtime_error {
        using std::runtime_error::runtime_error;
    };

    Fitter::Uncertainties_t photon_cb;
    Fitter::Uncertainties_t photon_taps;
    Fitter::Uncertainties_t proton_cb;
    Fitter::Uncertainties_t proton_taps;

    Constant();
    virtual ~Constant();

    Fitter::Uncertainties_t GetSigmas(const TParticle &particle) const override;

    /**
     * @brief Create a new instance and return a shared pointer to it
     * @return
     */
    static std::shared_ptr<Constant> make();
};

/**
 * @brief Simple kin fitter uncertainty model. has a fixed value for {cb,taps}{photon,proton}{theta,phi}, Energries are relative values and get multipied with the particle energy on GetSigmas()
 */
struct ConstantRelativeE : public Constant {
public:

    ConstantRelativeE();
    virtual ~ConstantRelativeE();

    Fitter::Uncertainties_t GetSigmas(const TParticle &particle) const override;

    /**
     * @brief Create a new, empty instance and return a shared pointer to it
     * @return
     */
    static std::shared_ptr<ConstantRelativeE> make();

    /**
     * @brief Create a new instance filled with global values determined from MC and return a shared ptr to it
     * @return
     */
    static std::shared_ptr<ConstantRelativeE> makeMCLongTarget();
};

struct ConstantRelativeEpow : public Constant {
public:
    double Eexp_cb   = 1.0;
    double Eexp_taps = 1.0;

    ConstantRelativeEpow();
    virtual ~ConstantRelativeEpow();

    Fitter::Uncertainties_t GetSigmas(const TParticle &particle) const override;

    /**
     * @brief Create a new, empty instance and return a shared pointer to it
     * @return
     */
    static std::shared_ptr<ConstantRelativeEpow> make();

    /**
     * @brief Create a new instance filled with global values determined from MC and return a shared ptr to it
     * @return
     */
    static std::shared_ptr<ConstantRelativeEpow> makeMCLongTarget();
};

/**
 * @brief Kin fitter uncertainties, uses histograms. Energy depenent values for each detector element. Histograms can be loaded from root files in setup database.
 */
class MCExtracted : public Fitter::UncertaintyModel {
public:
    struct angular_sigma {
        using Hist = std::shared_ptr<TH1D>;
        Hist p0 = nullptr;
        Hist p1 = nullptr;
        Hist p2 = nullptr;

        double GetSigma(const unsigned element, const double E) const;

        static double f(const double x, const double p0, const double p1, const double p2) noexcept;
        static double f_root(const double* x, const double* p) noexcept;

        static TF1* GetTF1(const std::string& name="SigmaFit");

        void Load(ant::WrapTFile& f, const std::string& prefix, const int bins);
        Hist LoadHist(ant::WrapTFile& f, const std::string& name, const int bins);

        angular_sigma();
        ~angular_sigma();

    };

protected:

    angular_sigma cb_sigma_theta;
    angular_sigma cb_sigma_phi;
    angular_sigma taps_sigma_theta;
    angular_sigma taps_sigma_phi;

    Fitter::Uncertainties_t GetSigmasProton(const TParticle &proton) const;
    Fitter::Uncertainties_t GetSigmasPhoton(const TParticle &photon) const;

public:

    struct Exception : std::runtime_error {
        using std::runtime_error::runtime_error;
    };

    MCExtracted();
    virtual ~MCExtracted();

    /**
     * @brief Load Sigmas from histograms in ROOT file
     * @param path Path to root file
     */
    void LoadSigmas(const std::string& path);

    Fitter::Uncertainties_t GetSigmas(const TParticle &particle) const override;

    /**
     * @brief Create an instance of this model and directly load sigmas from ROOT file of current setup
     * @return new instance
     */
    static std::shared_ptr<MCExtracted> makeAndLoad();

};

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
struct Optimized : Fitter::UncertaintyModel {

    double cb_photon_theta_const = 0.0;
    double cb_photon_theta_Sin   = 0.0;
    double cb_photon_phi         = 0.0;

    double cb_photon_E_rel = 0.0;
    double cb_photon_E_exp = 0.0;
    double cb_photon_E_lin = 0.0;

    Fitter::Uncertainties_t cb_proton = {};


    double taps_photon_E_rel = 0.0;
    double taps_photon_E_exp = 0.0;
    double taps_photon_E_lin = 0.0;
    double taps_photon_theta = 0.0;
    double taps_photon_phi   = 0.0;

    Fitter::Uncertainties_t taps_proton = {};

    struct Exception : std::runtime_error {
        using std::runtime_error::runtime_error;
    };

    Optimized();

    virtual ~Optimized();

    Fitter::Uncertainties_t GetSigmas(const TParticle& particle) const;
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
    const static std::string sepatator;

    void ReadToken(const std::string& token);
};

/**
 * @brief Values for the optimized model, obtainded from KinFitPi0 on data
 */
struct Optimized_Oli1 : Optimized {
    Optimized_Oli1();
};

}}}} // namespace ant::analysis::utils::UncertaintyModels