#pragma once

#include "analysis/physics/Physics.h"
#include "plot/PromptRandomHist.h"
#include "root-addons/cbtaps_display/TH2CB.h"
#include "base/WrapTTree.h"
#include "TTree.h"
#include "analysis/plot/PromptRandomHist.h"
#include "analysis/utils/TriggerSimulation.h"
#include <vector>

class TH1D;
class TTree;

namespace ant {
namespace analysis {
namespace physics {

class scratch_sobotzik_Pi0Calib : public Physics {
protected:

    struct MyTree : WrapTTree {


        ADD_BRANCH_T(double,E1)
        ADD_BRANCH_T(double,E2)
        ADD_BRANCH_T(double,M)
        ADD_BRANCH_T(double, Theta1)
        ADD_BRANCH_T(double, Theta2)
        ADD_BRANCH_T(double,Phi1)
        ADD_BRANCH_T(double,Phi2)
        ADD_BRANCH_T(double,OpeningAngle)
        ADD_BRANCH_T(double,ClusterSize1)
        ADD_BRANCH_T(double,ClusterSize2)
        ADD_BRANCH_T(double, ZVertex)
        ADD_BRANCH_T(int, ClusterNumber1)
        ADD_BRANCH_T(int, ClusterNumber2)
        ADD_BRANCH_T(double,true_E1)
        ADD_BRANCH_T(double,true_E2)
        ADD_BRANCH_T(double,true_openingangle)
        ADD_BRANCH_T(double, true_m)
        ADD_BRANCH_T(double, w)
    };



    TH1D* h_IM_All;

    TH3D* h_Meson_Energy_interval;
    TH3D* h_Meson_Energy_interval_30_Degree_Cut;

    TH2D* h_IM_CB_all;
    TH2D* h_IM_CB_Uncharged_No_Cut;
    TH2D* h_IM_CB_interval;
    TH2D* h_IM_CB_interval_Uncharged_No_Cut;
    TH2D* h_IM_CB_interval_Uncharged_30_Degree_Cut;
    TH2D* h_IM_CB_Angle_Energy;
    TH2D* h_IM_CB_Min_Opening_Angle;
    TH3D* h_IM_CB_Rec_vs_Gen_Opening_Angle;
    TH2D* h_IM_CB_Rec_vs_Gen_Opening_Angle_Deviation;
    TH2D* h_IM_CB_One_high_Photon;
    TH2D* h_IM_CB_AngleDeviation_Energy;
    TH3D* h_IM_CB_AngleDeviation_Photon_Meson_Energy;
//    TH2D* h_IM_CB_Rec_vs_Gen_Energie;
//    TH2D* h_IM_CB_Rec_Gen_Energie_Deviation;
    TH2D* h_IM_CB_Uncharged_30_Degree_Cut;
    //        TH3D* h_IM_CB_ZVertex;
    //        TH3D* h_IM_CB_ZVertex_interval;
    TH3D* h_IM_CB_ZVertex_interval_30_Degree_Cut;
    //        TH3D* h_IM_CB_Theta_Phi_Energy;
    TH3D* h_IM_CB_interval_Theta_Phi_Energy;
    TH1D* h_IM_CB_corr;
    TH1D* h_IM_TAPS;
    TH2D* h_IM_True_Opening_Angle;
    TH2D* h_IM_Rec_Opening_Angle;

    TH1D* h_Angle_CB;
    TH1D* h_Angle_TAPS;

    TH2D* h_ClusterHitTiming_CB;
    TH2D* h_ClusterHitTiming_TAPS;
    TH2D* h_IM_CB_ClusterSize3;
    TH2CB* h_cb;

    std::vector<TH2CB*> h_cbs_ClusterSize3;
    std::vector<TH2CB*> h_cbs_ClusterSize0;



    MyTree t;

    const interval<double> CaloEnergy_Window;
    ant::analysis::PromptRandom::Switch promptrandom;
    utils::TriggerSimulation triggersimu;

public:
    scratch_sobotzik_Pi0Calib(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};


}
}
}
