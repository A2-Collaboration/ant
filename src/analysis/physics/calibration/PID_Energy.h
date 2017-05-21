#pragma once

#include "analysis/physics/Physics.h"
#include "analysis/utils/fitter/KinFitter.h"
#include "analysis/utils/uncertainties/Optimized.h"
#include "analysis/utils/uncertainties/FitterSergey.h"
#include "analysis/plot/PromptRandomHist.h"
#include "analysis/utils/TriggerSimulation.h"

#include "TLorentzVector.h"

namespace ant {
namespace analysis {
namespace physics {

class PID_Energy : public Physics {

protected:
    TH2D* h_pedestals = nullptr;
    TH3D* h_bananas = nullptr;

    bool useHEP = false;
    // HEP related constants
    static constexpr bool PROBABILITY_CUT = true;
    static constexpr double PROBABILITY = .01;
    const unsigned MAX_GAMMA;
    // projection of the proton band
    static constexpr unsigned FIRST = 500;
    static constexpr unsigned LAST = 1000;

    struct PerChannel_t {
        TH2D* PedestalTiming = nullptr;
        TH1D* PedestalNoTiming = nullptr;
        TH2D* Banana = nullptr;
        TH2D* BananaRaw = nullptr;
        TH2D* BananaUnmatched = nullptr;
        TH1D* TDCMultiplicity;
        TH1D* QDCMultiplicity;
        PerChannel_t(analysis::HistogramFactory HistFac);
    };

    TH1D* h_BananaEntries;
    TH1D* h_BananaEntriesUnmatched;

    std::vector<PerChannel_t> h_perChannel;

    PromptRandom::Switch promptrandom;
    utils::TriggerSimulation triggersimu;
    utils::UncertaintyModelPtr model;
    utils::KinFitter kinfit;

    // containers for different multiplicities related to the HEP method
    TH2D* dEvE_all_combined;
    std::vector<TH2D*> dEvE_combined;
    TH2D* projections;
    std::vector<utils::KinFitter> kinfits;

    template<typename T>
    bool shift_right(std::vector<T>&);

    unsigned MinNGamma() const noexcept { return 2; }
    unsigned MaxNGamma() const noexcept { return MAX_GAMMA; }

public:

    PID_Energy(const std::string& name, OptionsPtr opts);

    static APLCON::Fit_Settings_t MakeFitSettings(unsigned);

    void ProcessHEP(const TEvent& event);
    bool find_best_comb(const TTaggerHit&, TCandidatePtrList&, TParticlePtr&);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
    virtual void Finish() override;
};

}}} // namespace ant::analysis::physics
