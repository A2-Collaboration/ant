#pragma once

#include "analysis/utils/Uncertainties.h"

class TH1D;
class TF1;

namespace ant {

class WrapTFile;

namespace analysis {
namespace utils {
namespace UncertaintyModels {

/**
 * @brief Kin fitter uncertainties, uses histograms. Energy dependent values for each detector element.
 * Histograms can be loaded from root files in setup database.
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

}}}} // namespace ant::analysis::utils::UncertaintyModels