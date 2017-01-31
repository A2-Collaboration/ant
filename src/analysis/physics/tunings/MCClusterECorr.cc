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
    const BinSettings bins_Ek(16*5,0,1600);

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

    if(event.Reconstructed().Candidates.size() != 1)
        return;

    auto& p_true = event.MCTrue().Particles.GetAll().front();
    auto& p_rec  = event.Reconstructed().Candidates.front();

    const auto cosThetaRec = cos(p_rec.Theta);
    const auto EkRec = p_rec.CaloEnergy;
    const auto EtrueErec = p_true->Ek()/EkRec;

    h_nFills_CB->Fill(EkRec, cosThetaRec, 1.0);
    h_EtrueErec_CB->Fill(EkRec, cosThetaRec, EtrueErec);

    h_nFills_TAPS->Fill(EkRec, cosThetaRec, 1.0);
    h_EtrueErec_TAPS->Fill(EkRec, cosThetaRec, EtrueErec);
}

void MCClusterECorr::Finish()
{
    h_EtrueErec_CB->Divide(h_nFills_CB);
    h_EtrueErec_TAPS->Divide(h_nFills_TAPS);
}

void MCClusterECorr::ShowResult()
{
    h_EtrueErec_CB->GetZaxis()->SetRangeUser(0.2,2);
    h_EtrueErec_TAPS->GetZaxis()->SetRangeUser(0.2,2);
    canvas(GetName()) << drawoption("colz")
                      << h_EtrueErec_CB
                      << h_nFills_CB
                      << h_EtrueErec_TAPS
                      << h_nFills_TAPS
                      << endc;
}

AUTO_REGISTER_PHYSICS(MCClusterECorr)