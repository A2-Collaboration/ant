#pragma once

#include "analysis/physics/Physics.h"
#include <map>
#include "base/WrapTTree.h"
#include "expconfig/detectors/Tagger.h"
#include "expconfig/detectors/PID.h"
#include "expconfig/detectors/TAPSVeto.h"
#include "expconfig/detectors/CB.h"
#include "analysis/utils/TriggerSimulation.h"
#include "analysis/plot/PromptRandomHist.h"
#include "analysis/utils/ProtonPhotonCombs.h"

#include "TLorentzVector.h"

class TH1D;

namespace ant {
namespace analysis {
namespace physics {

/**
 * @brief A class for analyzing the pi0->e+e-g decay
 *
 */
class Pi0Dalitz: public Physics {

protected:
    PromptRandom::Switch promptrandom;
    utils::TriggerSimulation triggersimu;
    using particle_comb_t = utils::ProtonPhotonCombs::comb_t;
    using particle_combs_t = utils::ProtonPhotonCombs::Combinations_t;

    std::shared_ptr<expconfig::detector::Tagger> tagger_detector;
    std::shared_ptr<expconfig::detector::CB> cb_detector;
    std::shared_ptr<expconfig::detector::PID> pid_detector;
    std::shared_ptr<expconfig::detector::TAPSVeto> veto_detector;

    static constexpr auto radtodeg = std_ext::radian_to_degree(1.0);
    enum parttype{en_p=0, en_ep, en_em, en_g};
    static const int nrPartTypes = 4;

    void CreateHistos();
    void DoTrueMCStuff(const std::vector<bool> &WhichMC, const std::vector<TParticlePtr> &trueparts);
    void DoMatchTrueRecoStuff(const TParticleList &allmcpart, const std::vector<TParticlePtr> &trueparts, const TCandidateList &recocands, std::vector<TParticlePtr> &matchrecopart);
    void DoTaggerStuff(const int cut, const TLorentzVector &g, const double &time, const double &tw);
    void DoTriggerStuff(const int cut, const double &tw);
    void DoRecoCandStuff(const int cut, const TCandidateList &recocands, particle_combs_t ppcomb, const std::vector<TParticlePtr> &recmatparts, const std::vector<bool> &WhichMC, const TLorentzVector &ig, const double &tw);

    //-- Histograms
    static const int nrSel = 5;
    //--- MCtrue
    TH1D *h_IMeegTrue, *h_IMggTrue, *h_RecoTrueMatch;
    TH2D *h_RecoTrueAngle, *h_ThetavsEnergy_MCTrue[nrPartTypes];
    //--- Trigger
    TH1D *h_CBEsum[nrSel];
    //--- Tagger
    TH1D *h_TaggTimeuw[nrSel], *h_TaggTimeww[nrSel],*h_TaggPhEnww[nrSel], *h_TaggPhEnuw[nrSel];
    //--- Reconstructed candidate info
    TH2D *h_PIDMultUsed, *h_VetoMultUsed;
    TH1D *h_IMeegReco[nrSel], *h_IMggReco[nrSel], *h_MMpReco[nrSel], *h_OpAngpphReco[nrSel];
    TH2D *h_PIDEvsT[nrSel], *h_TVetoEvsT[nrSel], *h_CBEvsT[nrSel], *h_TAPSEvsT[nrSel], *h_CBEvsNrCr[nrSel], *h_TAPSEvsNrCr[nrSel];
    TH2D *h_TimeCBvsPID[nrSel], *h_EnergyCBvsPID[nrSel], *h_TimeTAPSvsTVeto[nrSel], *h_EnergyTAPSvsTVeto[nrSel];
    TH2D *h_ThPhCBvsPID[nrSel], *h_ThPhTAPSvsTVeto[nrSel], *h_ThetavsEnergy[nrSel];
    TH2D *h_EnergyCBvsPID_RecMat[nrSel][nrPartTypes], *h_EnergyTAPSvsTVeto_RecMat[nrSel][nrPartTypes];
    TH2D *h_ThetavsEnergy_RecMat[nrSel][nrPartTypes], *h_CBEvsNrCr_RecMat[nrSel][nrPartTypes], *h_TAPSEvsNrCr_RecMat[nrSel][nrPartTypes];
    TH1D *h_OpAngpphReco_RecMat[nrSel];

public:
    Pi0Dalitz(const std::string& name, OptionsPtr opts);
    virtual ~Pi0Dalitz() {}

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override;
    virtual void ShowResult() override;

};

}
}
}
