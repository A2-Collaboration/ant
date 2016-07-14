#pragma once

#include "tree/TParticle.h"
#include "base/printable.h"

#include "TH1D.h"

class TRandom;

namespace ant {

class WrapTFile;
class Interpolator2D;

namespace analysis {
namespace utils {

/**
 * @brief Uncertainties for E, theta, and phi
 */
struct Uncertainties_t {
    double sigmaE     = {};
    double sigmaTheta = {};
    double sigmaPhi   = {};

    Uncertainties_t() = default;
    Uncertainties_t(const double E, const double Theta, const double Phi) : sigmaE(E), sigmaTheta(Theta), sigmaPhi(Phi) {}

    bool operator==(const Uncertainties_t& other) const noexcept {
        return sigmaE == other.sigmaE && sigmaTheta == other.sigmaTheta && sigmaPhi == other.sigmaPhi;
    }
};

/**
 * @brief Virtual base class for different Uncertainty Models for kion fitter.
 *        Derive and implement the GetSigmas() method.
 * @see UncertaintyModels::Constant
 * @see UncertaintyModels::MCExtracted
 */
class UncertaintyModel {
public:
    virtual ~UncertaintyModel();
    virtual Uncertainties_t GetSigmas(const TParticle& particle) const =0;
};
using UncertaintyModelPtr = std::shared_ptr<const UncertaintyModel>;

namespace UncertaintyModels {

/**
 * @brief Constant/static kin fitter uncertainty model. has a fixed value for {cb,taps}{photon,proton}{E,theta,phi}
 */
struct Constant : public UncertaintyModel {
public:

    struct Exception : std::runtime_error {
        using std::runtime_error::runtime_error;
    };

    Uncertainties_t photon_cb;
    Uncertainties_t photon_taps;
    Uncertainties_t proton_cb;
    Uncertainties_t proton_taps;

    Constant();
    virtual ~Constant();

    Uncertainties_t GetSigmas(const TParticle &particle) const override;

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

    Uncertainties_t GetSigmas(const TParticle &particle) const override;

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

    Uncertainties_t GetSigmas(const TParticle &particle) const override;

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
class MCExtracted : public UncertaintyModel {
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

    Uncertainties_t GetSigmasProton(const TParticle &proton) const;
    Uncertainties_t GetSigmasPhoton(const TParticle &photon) const;

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

    Uncertainties_t GetSigmas(const TParticle &particle) const override;

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

    struct Exception : std::runtime_error {
        using std::runtime_error::runtime_error;
    };

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
    Optimized_Oli1();
};

/**
 * @brief Uncertainties from Patrik Adlarson for MC Smearing
 */
struct MCSmearingAdlarson : public UncertaintyModel {
public:

    MCSmearingAdlarson();
    virtual ~MCSmearingAdlarson();

    Uncertainties_t GetSigmas(const TParticle &particle) const override;

    /**
     * @brief Create a new instance and return a shared pointer to it
     * @return
     */
    static std::shared_ptr<MCSmearingAdlarson> make();

protected:
    std::unique_ptr<TRandom> rng;
};

/**
 * @brief Uncertainties with interpolated surfaces in (E,theta) plane,
 * determined with iterative procedure
 */
struct Interpolated : public UncertaintyModel, public ant::printable_traits {
public:

    Interpolated(UncertaintyModelPtr starting_uncertainty_);
    virtual ~Interpolated();

    Uncertainties_t GetSigmas(const TParticle &particle) const override;

    void LoadSigmas(const std::string& filename);

    struct Exception : std::runtime_error {
        using std::runtime_error::runtime_error;
    };

    static std::shared_ptr<Interpolated> makeAndLoad(UncertaintyModelPtr default_model);

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
    UncertaintyModelPtr starting_uncertainty;

    bool loaded_sigmas = false;

    struct EkThetaPhi : ant::printable_traits {

        ClippedInterpolatorWrapper E;
        ClippedInterpolatorWrapper Theta;
        ClippedInterpolatorWrapper Phi;

        Uncertainties_t GetUncertainties(const TParticle& particle) const;

        void Load(ant::WrapTFile& file, const std::string& prefix);
        static std::unique_ptr<const Interpolator2D> LoadInterpolator(ant::WrapTFile& file, const std::string& prefix);

        std::ostream& Print(std::ostream& stream) const override;
    };

    EkThetaPhi cb_photon;
    EkThetaPhi taps_photon;
    EkThetaPhi cb_proton;
    EkThetaPhi taps_proton;
};

}}}} // namespace ant::analysis::utils::UncertaintyModels
