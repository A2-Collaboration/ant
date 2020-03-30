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

    struct RecoCandInfo
    {
        double CaloE; double VetoE;
        double Theta; double Phi; double CaloT;
        int CaloCluSize; bool InCB;
        double VetoTheta; double VetoPhi; double VetoT;
        int VetoEl; int VetoElClean;
        bool InPID;

        void PrintInfo(){
            printf("Ecalo=%.2f, Eveto=%.2f, CalTheta=%.2f, CalPhi=%.2f, CalT=%.2f, CluSize=%d, In CB=%d\n", CaloE, VetoE, Theta*radtodeg, Phi*radtodeg, CaloT, CaloCluSize, InCB);
            printf("VetoTheta=%.2f, VetoPhi=%.2f, VetoT=%.2f, VetoElNr=%d, VetoClElNr=%d, In PID=%d\n\n", VetoTheta*radtodeg, VetoPhi*radtodeg, VetoT, VetoEl, VetoElClean, InPID);
        }
    };

    void CreateHistos();
    void CheckVetoUsage(const int elnr, const bool inpid, int &elnrcleaned, std::map<int,int> &pidfreq, std::map<int,int> &vetofreq);
    void FillRecoInfo(RecoCandInfo *rc, const std::vector<double> doublelist, const std::vector<int> intlist, const std::vector<bool> boollist);
    void DoTrueMCStuff(const int WhichMC, const TLorentzVector &ep, const TLorentzVector &em, const std::vector<TLorentzVector> &g);
    void DoTaggerStuff(const int cut, const TLorentzVector &g, const double &time, const double &tw);
    void DoTriggerStuff(const int cut, const double &tw);
    void DoRecoCandStuff(const int cut, const std::vector<RecoCandInfo> &recocandinfo, std::map<int,int> &pidfreq, std::map<int,int> &vetofreq, particle_combs_t ppcomb, const TLorentzVector &ig, const double &tw);

    // histograms
    TH1D *h_IMeegTrue, *h_IsPi0eegMC;
    TH2D *h_PIDMultUsed, *h_VetoMultUsed;
    static const int nrSel = 5;
    TH1D *h_CBEsum[nrSel];
    TH1D *h_TaggTimeuw[nrSel], *h_TaggTimeww[nrSel],*h_TaggPhEnww[nrSel], *h_TaggPhEnuw[nrSel];
    TH1D *h_IMeegReco[nrSel], *h_IMggReco[nrSel], *h_MMpReco[nrSel], *h_OpAngpphReco[nrSel];
    TH2D *h_PIDEvsT[nrSel], *h_TVetoEvsT[nrSel], *h_CBEvsT[nrSel], *h_TAPSEvsT[nrSel], *h_CBEvsNrCr[nrSel], *h_TAPSEvsNrCr[nrSel];
    TH2D *h_TimeCBvsPID[nrSel], *h_EnergyCBvsPID[nrSel], *h_TimeTAPSvsTVeto[nrSel], *h_EnergyTAPSvsTVeto[nrSel];
    TH2D *h_ThPhCBvsPID[nrSel], *h_ThPhTAPSvsTVeto[nrSel];

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
