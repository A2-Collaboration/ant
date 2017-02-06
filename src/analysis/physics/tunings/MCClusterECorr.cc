#include "MCClusterECorr.h"

#include "base/Logger.h"
#include "expconfig/ExpConfig.h"

#include "plot/root_draw.h"



// use some namespaces (remember: only in implementation (aka .cc) files
// those statements are recommended to keep the following code not so namespace-clobbered
using namespace std;
using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::physics;


MCClusterECorr::MCClusterECorr(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    CB(Detector_t::Type_t::CB,     HistFac,  {35, cos(degree_to_radian(160.0)), cos(degree_to_radian(20.0))}),
    TAPS(Detector_t::Type_t::TAPS, HistFac,  {10, cos(degree_to_radian( 20.0)), cos(degree_to_radian( 0.0))})
{
    const BinSettings bins_Ek(16*4,0,1600);

    h_nCaloClusters = HistFac.makeTH2D("nCaloClusters","E_{kin}^{rec} / MeV","nCaloClusters",
                                       bins_Ek, BinSettings(5),
                                       "h_nCaloClusters");

    h_LostMCTrue     = HistFac.makeTH2D("nCaloClusters_{CB}!=1 && nCaloClusters_{TAPS}!=1","E_{kin}^{true} / MeV","cos #theta^{true}",
                                        bins_Ek, BinSettings(40,-1, 1), "h_LostMCTrue");
}

MCClusterECorr::CBTAPS_t::CBTAPS_t(Detector_t::Type_t type,
                                   const HistogramFactory& histFac,
                                   const BinSettings& bins_cosTheta) :
    Detector(dynamic_pointer_cast<ClusterDetector_t>(ExpConfig::Setup::GetDetector(type))),
    HistFac(Detector_t::ToString(type), histFac, Detector_t::ToString(type))
{
    const BinSettings bins_Ek(16*4,0,1600);

    h_nFills      = HistFac.makeTH2D("nFills","E_{kin}^{rec} / MeV","cos #theta^{rec}",
                                     bins_Ek, bins_cosTheta, "h_nFills");
    h_EtrueErec   = HistFac.makeTH2D("<E^{true}/E^{rec}>",  "E_{kin}^{rec} / MeV","cos #theta^{rec}",
                                     bins_Ek, bins_cosTheta, "h_EtrueErec");
    h_EtrueErec->SetBit(TH1::kIsAverage); // pretty important for Ant-hadd

    h_EtrueErec_3D   = HistFac.makeTH3D("<E^{true}/E^{rec}>",  "E_{kin}^{rec} / MeV","cos #theta^{rec}","Etrue/Erec",
                                        bins_Ek, bins_cosTheta, BinSettings(40,0.5,2), "h_EtrueErec_3D");

    h_ErecEtrue_elements   = HistFac.makeTH2D("E^{rec}/E^{true}",  "E^{rec}/E^{true}","Element",
                                              BinSettings(100,0,2),
                                              BinSettings(Detector->GetNChannels()), "h_ErecEtrue_elements");
}

void MCClusterECorr::CBTAPS_t::Fill(const TCluster& caloCluster, double Etrue) const
{
    if(caloCluster.DetectorType != Detector->Type)
        return;

    if(Detector->GetClusterElement(caloCluster.CentralElement)->TouchesHole)
        return;

    const auto cosThetaRec = cos(caloCluster.Position.Theta());
    const auto EkRec = caloCluster.Energy;
    const auto EtrueErec = Etrue/EkRec;

    h_nFills->Fill(EkRec, cosThetaRec, 1.0);
    h_EtrueErec->Fill(EkRec, cosThetaRec, EtrueErec);
    h_EtrueErec_3D->Fill(EkRec, cosThetaRec, EtrueErec);
    h_ErecEtrue_elements->Fill(EkRec/Etrue, caloCluster.CentralElement);
}

void MCClusterECorr::CBTAPS_t::Finish() const
{
    auto divide_hist = [] (TH2* dst, const TH2* src) {
        for(int binx=0; binx<=src->GetNbinsX()+1; ++binx) {
            for(int biny=0; biny<=src->GetNbinsY()+1; ++biny) {
                const auto divisor = src->GetBinContent(binx,biny);
                const auto res = dst->GetBinContent(binx,biny)/divisor;
                dst->SetBinContent(binx, biny, isfinite(res) ? res : 1.0);
            }
        }
    };

    divide_hist(h_EtrueErec, h_nFills);
}

void ant::analysis::physics::MCClusterECorr::CBTAPS_t::Draw(canvas& c) const
{
    h_EtrueErec->GetZaxis()->SetRangeUser(0.2,2);
    c << h_EtrueErec
      << h_nFills
      << h_ErecEtrue_elements;
}

void MCClusterECorr::ProcessEvent(const TEvent& event, manager_t&)
{
    // we can only run on single particle gun MC data
    // preferably photons...
    if(event.MCTrue().Particles.GetAll().size() != 1)
        return;
    auto& p_true = event.MCTrue().Particles.GetAll().front();

    // determine number of calo clusters
    auto nCaloClusters_CB = 0;
    auto nCaloClusters_TAPS = 0;
    for(auto& cl : event.Reconstructed().Clusters) {
        if(cl.DetectorType & Detector_t::Type_t::CB)
            nCaloClusters_CB++;
        else if(cl.DetectorType & Detector_t::Type_t::TAPS)
            nCaloClusters_TAPS++;
    }

    // fill some check histogram
    for(auto& cl : event.Reconstructed().Clusters) {
        if(cl.DetectorType & Detector_t::Any_t::Calo)
            h_nCaloClusters->Fill(cl.Energy, nCaloClusters_CB + nCaloClusters_TAPS);
    }

    // only consider events with exactly one cluster reconstructed
    if((nCaloClusters_CB != 1) && (nCaloClusters_TAPS != 1)) {
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

    if(nCaloClusters_CB == 1)
        CB.Fill(*caloCluster, p_true->Ek());

    if(nCaloClusters_TAPS ==1)
        TAPS.Fill(*caloCluster, p_true->Ek());
}

void MCClusterECorr::Finish()
{
    CB.Finish();
    TAPS.Finish();
}

namespace ant {
template<typename canvas>
canvas& operator<<(canvas&& c, const MCClusterECorr::CBTAPS_t& cbtaps) {
    cbtaps.Draw(c);
    return c;
}
}

void MCClusterECorr::ShowResult()
{
    canvas(GetName())
            << drawoption("colz")
            << padoption::LogZ << h_nCaloClusters
            << h_LostMCTrue
            << endr
            << CB
            << endr
            << TAPS
            << endc;
}





AUTO_REGISTER_PHYSICS(MCClusterECorr)


