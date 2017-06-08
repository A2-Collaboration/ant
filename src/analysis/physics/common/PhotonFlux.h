#pragma once

#include "analysis/physics/Physics.h"
#include "base/piecewise_interval.h"
#include "base/WrapTTree.h"

class TH1D;

namespace ant {
namespace analysis {
namespace physics {


class PhotonFlux: public Physics {
protected:
    const std::shared_ptr<TaggerDetector_t> tagger;

    unsigned seenScalerBlocks = 0;
    unsigned nchannels = 0;

    double time = 0;

    const double targetDensity =  4.249E-7;  // TargetDensity in microbarn^-1


public:
    TH1D* ScalerCounts = nullptr;
    TH1D* TaggEff      = nullptr;

    TH1D* Flux         = nullptr;
    TH1D* FluxLTcor    = nullptr;
    TH1D* lifetime     = nullptr;

    TH1D* IntLumi      = nullptr;
    TH1D* IntLumiLTcor = nullptr;
    TH1D* Lumi         = nullptr;

    TH1D* info         = nullptr;

    PhotonFlux(const std::string& name, OptionsPtr opts=nullptr);
    virtual ~PhotonFlux();

    virtual void ProcessEvent(const TEvent& ev, manager_t&) override;
    virtual void Finish() override;
    virtual void ShowResult() override;

    void processBlock();
};

}}} // namespace ant::analysis::physics
