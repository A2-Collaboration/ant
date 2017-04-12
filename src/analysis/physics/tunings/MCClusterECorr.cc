#include "MCClusterECorr.h"

#include "base/Logger.h"
#include "expconfig/ExpConfig.h"

#include "utils/ParticleTools.h"



// use some namespaces (remember: only in implementation (aka .cc) files
// those statements are recommended to keep the following code not so namespace-clobbered
using namespace std;
using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::physics;

const BinSettings MCClusterECorr::bins_Ek     = BinSettings(160,0,1600);
const BinSettings MCClusterECorr::bins_clSize = BinSettings(30);

MCClusterECorr::MCClusterECorr(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    CBThetaWindow(degree_to_radian(50.0), degree_to_radian(180.0-50.0)),
    CBHemisphereGap({degree_to_radian(interval<double>::CenterWidth(0.0,40.0)),degree_to_radian(interval<double>::CenterWidth(180.0,40.0))}),
    CB(Detector_t::Type_t::CB,     HistFac),
    TAPS(Detector_t::Type_t::TAPS, HistFac)
{

    h_nCaloClusters = HistFac.makeTH2D("nCaloClusters","E_{kin}^{rec} / MeV","nCaloClusters",
                                       bins_Ek, BinSettings(5),
                                       "h_nCaloClusters");

    h_LostMCTrue     = HistFac.makeTH2D("nCaloClusters_{CB}!=1 && nCaloClusters_{TAPS}!=1","E_{kin}^{true} / MeV","cos #theta^{true}",
                                        bins_Ek, BinSettings(40,-1, 1), "h_LostMCTrue");
}

MCClusterECorr::CBTAPS_t::CBTAPS_t(Detector_t::Type_t type,
                                   const HistogramFactory& histFac) :
    Type(type),
    HistFac(Detector_t::ToString(Type), histFac, Detector_t::ToString(Type))
{
    auto det = ExpConfig::Setup::GetDetector(Type);

    h_nFills      = HistFac.makeTH2D("nFills","E_{kin}^{rec} / MeV","Cluster Size",
                                     bins_Ek, bins_clSize, "h_nFills");
    h_EtrueErec   = HistFac.makeTH2D("<E^{true}/E^{rec}>",  "E_{kin}^{rec} / MeV","Cluster Size",
                                     bins_Ek, bins_clSize, "h_EtrueErec");
    h_EtrueErec->SetBit(TH1::kIsAverage); // pretty important for Ant-hadd

    h_EtrueErec_3D   = HistFac.makeTH3D("<E^{true}/E^{rec}>",  "E_{kin}^{rec} / MeV","Cluster Size","Etrue/Erec",
                                        bins_Ek, bins_clSize, BinSettings(40,0.5,2), "h_EtrueErec_3D");

    h_ErecEtrue_elements   = HistFac.makeTH2D("E^{rec}/E^{true}",  "E^{rec}/E^{true}","Element",
                                              BinSettings(100,0,2),
                                              BinSettings(det->GetNChannels()), "h_ErecEtrue_elements");
}

void MCClusterECorr::CBTAPS_t::Fill(const TCluster& caloCluster, double Etrue) const
{
    if(caloCluster.DetectorType != Type)
        return;

    if(caloCluster.HasFlag(TCluster::Flags_t::TouchesHoleCentral))
        return;

    const auto EkRec = caloCluster.Energy;
    const auto EtrueErec = Etrue/EkRec;
    const auto clSise = caloCluster.Hits.size();

    h_nFills->Fill(EkRec, clSise, 1.0);
    h_EtrueErec->Fill(EkRec, clSise, EtrueErec);
    h_EtrueErec_3D->Fill(EkRec, clSise, EtrueErec);
    h_ErecEtrue_elements->Fill(EkRec/Etrue, caloCluster.CentralElement);
}

void MCClusterECorr::CBTAPS_t::Finish() const
{
    auto divide_hist = [] (TH2* dst, const TH2* src, const double min=0.0) {
        for(int binx=0; binx<=src->GetNbinsX()+1; ++binx) {
            for(int biny=0; biny<=src->GetNbinsY()+1; ++biny) {
                const auto divisor = src->GetBinContent(binx,biny);
                const auto res = dst->GetBinContent(binx,biny)/divisor;
                dst->SetBinContent(binx, biny, isfinite(res) && (divisor > min) ?  res : -1.0);
            }
        }
    };

    divide_hist(h_EtrueErec, h_nFills, 0.0);
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
    auto mctrue_particles = utils::ParticleTypeList::Make(event.MCTrue().ParticleTree);

    // we can only run on single particle gun MC data
    // preferably photons...
    if(mctrue_particles.GetAll().size() != 1)
        return;
    auto& p_true = mctrue_particles.GetAll().front();

    // determine number of calo clusters
    auto nCaloClusters_CB = 0;
    auto nCaloClusters_TAPS = 0;
    for(auto& cl : event.Reconstructed().Clusters) {
        if(cl.DetectorType == Detector_t::Type_t::CB)
            nCaloClusters_CB++;
        else if(cl.DetectorType == Detector_t::Type_t::TAPS)
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

    if(nCaloClusters_CB == 1) {
        if(CBThetaWindow.Contains(caloCluster->Position.Theta())
                && !CBHemisphereGap.Contains(caloCluster->Position.Phi()))
            CB.Fill(*caloCluster, p_true->Ek());
    }

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


