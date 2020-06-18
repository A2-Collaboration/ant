#pragma once

#include "physics/Physics.h"
#include <map>
#include "base/WrapTTree.h"
#include "utils/TriggerSimulation.h"
#include "analysis/plot/PromptRandomHist.h"
#include "expconfig/detectors/Tagger.h"
#include "TLorentzVector.h"


class TH1D;
using namespace std;

namespace ant {
namespace analysis {
namespace physics {

/**
 * @brief A gp->pw->pi0g class
 *
 */
class scratch_lheijken_gpwppi0g: public Physics {
protected:

    PromptRandom::Switch promptrandom;
    utils::TriggerSimulation triggersimu;
    std::shared_ptr<expconfig::detector::Tagger> tagger_detector;
    unsigned int nTagger;

    struct str_wpi0gO3g
    {
        vector<TLorentzVector> LV3g;
        TLorentzVector LVw;
        TLorentzVector LVp;
        double avgT3g;
        bool onlyinCB;
    };

    struct str_recocand
    {
        TCandidatePtrList neutral;
        TCandidatePtrList charged;
        vector<double> NeuCanCaloE;
        vector<double> NeuCanTheta;
        vector<double> NeuCanPhi;
        vector<double> NeuCanTime;
        vector<double> NeuCanCluSize;
        vector<bool> NeuCanInCB;
        vector<double> ChaCanCaloE;
        vector<double> ChaCanVetoE;
        vector<double> ChaCanTheta;
        vector<double> ChaCanPhi;
        vector<double> ChaCanTime;
        vector<double> ChaCanCluSize;
        vector<bool> ChaCanInCB;
        vector<bool> ChaCanInPID;
    };

    //--- MC True
    TH1D *hTrueGammaE, *hTrueIMggg, *hTrueMMp;
    TH2D *hTrueThevsEg, *hTrueThevsPhig, *hTrueThevsEw, *hTrueThevsPhiw;
    //--- Checks
    TH1D *hPromRandWei;
    TH2D *hTaggTimeChannel, *hTaggCorTimeChannel;
    //--- All candidate info
    TH2D *hNeuCanThevsCaloE, *hNeuCanThevsPhi;
    TH1D *hNeuCanCBTime, *hNeuCanTAPSTime, *hNeuCanCBCluSize, *hNeuCanTAPSCluSize;
    TH2D *hChaCanThevsCaloE,*hChaCanThevsVetoE, *hChaCanThevsPhi;
    TH1D *hChaCanCBTime, *hChaCanTAPSTime, *hChaCanCBCluSize, *hChaCanTAPSCluSize;
    //--- w->pi0g stuff
    TH1D *hO3gOCB_time,*hO3gCBTA_time, *hO3gOCB_IM, *hO3gCBTA_IM, *hO3gOCB_MM, *hO3gCBTA_MM;
    TH2D *hO3g_gThevsE, *hO3g_gThevsPhi, *hO3g_wThevsE, *hO3g_wThevsPhi;
    TH2D *hO3gOCB_IMgggVsTagCh, *hO3gCBTA_IMgggVsTagCh;
    TH2D *hO3gOCB_TagTimeVsChan, *hO3gCBTA_TagTimeVsChan, *hO3gOCB_TagCorTimeVsChan, *hO3gCBTA_TagCorTimeVsChan;
    TH1D *hO3gO1p_IMggg, *hO3gO1p_MM;
    TH2D *hO3gO1p_IMggg_Thggg;

protected:
    static constexpr auto radtodeg = std_ext::radian_to_degree(1.0);
    void CreateHistos();
    void FillO3gHistos(const str_wpi0gO3g strwpi0g, const TLorentzVector initphoton, const TTaggerHit th, const double ctt, const double tw);
    void FillRecoCandHistos(const str_recocand rc, const double tw);

public:
    scratch_lheijken_gpwppi0g(const std::string& name, OptionsPtr opts);
    virtual ~scratch_lheijken_gpwppi0g() {}

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override;
    virtual void ShowResult() override;

};

}
}
}
