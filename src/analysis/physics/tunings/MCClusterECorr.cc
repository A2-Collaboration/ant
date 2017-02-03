#include "MCClusterECorr.h"

#include "base/Logger.h"

#include "plot/root_draw.h"



// use some namespaces (remember: only in implementation (aka .cc) files
// those statements are recommended to keep the following code not so namespace-clobbered
using namespace std;
using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::physics;


MCClusterECorr::MCClusterECorr(const string& name, OptionsPtr opts) :
    Physics(name, opts)
{

    // copied/synced with TwoPi0_MCSmearing (but finer Ek binning)
    const BinSettings bins_cosTheta_CB  (35, cos(degree_to_radian(160.0)), cos(degree_to_radian(20.0)));
    const BinSettings bins_cosTheta_TAPS(10, cos(degree_to_radian( 20.0)), cos(degree_to_radian( 0.0)));
    const BinSettings bins_Ek(16*2,0,1600);

    h_nCaloClusters = HistFac.makeTH2D("nCaloClusters","E_{kin}^{rec} / MeV","nCaloClusters",
                                       bins_Ek,BinSettings(5),
                                       "h_nCaloClusters");

    h_LostMCTrue     = HistFac.makeTH2D("nCaloClusters!=1","E_{kin}^{true} / MeV","cos #theta^{true}",
                                        bins_Ek, bins_cosTheta_CB, "h_LostMCTrue");


    h_nFills_CB      = HistFac.makeTH2D("nFills CB","E_{kin}^{rec} / MeV","cos #theta^{rec}",
                                        bins_Ek, bins_cosTheta_CB, "h_nFills_CB");
    h_EtrueErec_CB   = HistFac.makeTH2D("<E^{true}/E^{rec}> CB",  "E_{kin}^{rec} / MeV","cos #theta^{rec}",
                                        bins_Ek, bins_cosTheta_CB, "h_EtrueErec_CB");
    h_EtrueErec_CB->SetBit(TH1::kIsAverage);

    h_nFills_TAPS    = HistFac.makeTH2D("nFills TAPS","E_{kin}^{rec} / MeV","cos #theta^{rec}",
                                        bins_Ek, bins_cosTheta_TAPS, "h_nFills_TAPS");
    h_EtrueErec_TAPS = HistFac.makeTH2D("<E^{true}/E^{rec}> TAPS","E_{kin}^{rec} / MeV","cos #theta^{rec}",
                                        bins_Ek, bins_cosTheta_TAPS, "h_EtrueErec_TAPS");
    h_EtrueErec_TAPS->SetBit(TH1::kIsAverage);
}

void MCClusterECorr::ProcessEvent(const TEvent& event, manager_t&)
{
    // we can only run on single particle gun MC data
    // preferably photons...
    if(event.MCTrue().Particles.GetAll().size() != 1)
        return;
    auto& p_true = event.MCTrue().Particles.GetAll().front();

    // determine number of calo clusters
    auto nCaloClusters = 0;
    for(auto& cl : event.Reconstructed().Clusters) {
        if(cl.DetectorType & Detector_t::Any_t::Calo)
            nCaloClusters++;
    }

    // fill some check histogram
    for(auto& cl : event.Reconstructed().Clusters) {
        if(cl.DetectorType & Detector_t::Any_t::Calo)
            h_nCaloClusters->Fill(cl.Energy, nCaloClusters);
    }

    // only consider events with exactly one cluster reconstructed
    if(nCaloClusters != 1) {
        // fill check histogram
        h_LostMCTrue->Fill(p_true->Ek(), cos(p_true->Theta()), 1.0);
        return;
    }

    // find that one CaloCluster
    TClusterPtr caloCluster;
    for(auto& cl : event.Reconstructed().Clusters.get_iter()) {
        if(cl->DetectorType & Detector_t::Any_t::Calo)
            caloCluster = cl;
    }

    const auto cosThetaRec = cos(caloCluster->Position.Theta());
    const auto EkRec = caloCluster->Energy;
    const auto EtrueErec = p_true->Ek()/EkRec;

    h_nFills_CB->Fill(EkRec, cosThetaRec, 1.0);
    h_EtrueErec_CB->Fill(EkRec, cosThetaRec, EtrueErec);

    h_nFills_TAPS->Fill(EkRec, cosThetaRec, 1.0);
    h_EtrueErec_TAPS->Fill(EkRec, cosThetaRec, EtrueErec);
}

void MCClusterECorr::Finish()
{
    auto divide_hist = [] (TH2* dst, const TH2* src) {
        for(int binx=0; binx<=src->GetNbinsX(); ++binx) {
            for(int biny=0; biny<=src->GetNbinsY(); ++biny) {
                const auto divisor = src->GetBinContent(binx,biny);
                const auto res = dst->GetBinContent(binx,biny)/divisor;
                dst->SetBinContent(binx, biny, isfinite(res) ? res : 1.0);
            }
        }
    };

    divide_hist(h_EtrueErec_CB,   h_nFills_CB);
    divide_hist(h_EtrueErec_TAPS, h_nFills_TAPS);
}

void MCClusterECorr::ShowResult()
{
    h_EtrueErec_CB->GetZaxis()->SetRangeUser(0.2,2);
    h_EtrueErec_TAPS->GetZaxis()->SetRangeUser(0.2,2);
    canvas(GetName()) << drawoption("colz")
                      << padoption::LogZ << h_nCaloClusters
                      << h_LostMCTrue
                      << endr
                      << h_EtrueErec_CB
                      << h_nFills_CB
                      << endr
                      << h_EtrueErec_TAPS
                      << h_nFills_TAPS
                      << endc;
}

AUTO_REGISTER_PHYSICS(MCClusterECorr)