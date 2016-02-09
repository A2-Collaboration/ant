#pragma once

#include "analysis/physics/Physics.h"
#include "base/interval.h"
#include "plot/PromptRandomHist.h"
#include "utils/KinFitter.h"
#include "expconfig/detectors/TAPS.h"

class TH1D;
class TH2D;

namespace ant {
namespace analysis {
namespace physics {

class EtapProton : public Physics {

protected:
    PiecewiseInterval<unsigned> multiplicities;

    TH1D* steps;

    TTree* tree = nullptr;

    unsigned  b_nCB = 0;
    unsigned  b_nTAPS = 0;
    double    b_CBAvgTime  = 0.0;
    double    b_CBSumVetoE  = 0.0;

    TLorentzVector b_PhotonSum;
    double b_ProtonCopl;
    double b_ProtonBeta;
    double b_ProtonToF;

    PromptRandom::Switch promptrandom;

    std::vector<std::unique_ptr<utils::KinFitter>> fitters;
    double b_FitChi2;
    unsigned b_FitStatus;
    unsigned b_NFitIterations;
    double b_TaggW;
    double b_TaggE;
    double b_TaggT;
    unsigned b_TaggCh;
    TLorentzVector b_Missing;


    double b_FittedTaggE;
    TLorentzVector b_FittedProton;
    std::vector<TLorentzVector> b_FittedPhotons;
    TLorentzVector b_FittedPhotonSum;
    double b_FittedProtonCopl;

    std::shared_ptr<expconfig::detector::TAPS> taps_detector;

public:

    EtapProton(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};

}
}
}
