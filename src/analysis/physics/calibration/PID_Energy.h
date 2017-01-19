#pragma once

#include "analysis/physics/Physics.h"
#include "analysis/utils/Fitter.h"
#include "analysis/utils/uncertainties/Optimized.h"
#include "analysis/utils/uncertainties/FitterSergey.h"
#include "analysis/plot/PromptRandomHist.h"

#include "TLorentzVector.h"

namespace ant {
namespace analysis {
namespace physics {

class PID_Energy : public Physics {

protected:
    TH2D* h_pedestals = nullptr;
    TH3D* h_bananas = nullptr;
    TH2D* h_mip = nullptr;

    bool useMIP = false;

    struct PerChannel_t {
        TH2D* PedestalTiming = nullptr;
        TH1D* PedestalNoTiming = nullptr;
        TH2D* Banana = nullptr;
        TH2D* BananaRaw = nullptr;
        TH1D* TDCMultiplicity;
        TH1D* QDCMultiplicity;
        PerChannel_t(analysis::HistogramFactory HistFac);
    };

    std::vector<PerChannel_t> h_perChannel;

    PromptRandom::Switch promptrandom;
    utils::UncertaintyModelPtr model;
    utils::KinFitter kinfit;

    template<typename T>
    void shift_right(std::vector<T>&);

public:

    PID_Energy(const std::string& name, OptionsPtr opts);

    static APLCON::Fit_Settings_t MakeFitSettings(unsigned);

    bool doFit_checkProb(const TTaggerHit& taggerhit,
                         const TParticlePtr proton,
                         const TParticleList photons,
                         double& best_prob_fit,
                         TParticleList& fit_photons);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};

}}} // namespace ant::analysis::physics
