#include "MCEnergyThresholds.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

#include "MCEnergyThresholds.h"

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

const BinSettings MCEnergyThresholds::bins_Ek     = BinSettings(160,0,1600);
const BinSettings MCEnergyThresholds::bins_clSize = BinSettings(30);

MCEnergyThresholds::MCEnergyThresholds(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    CBThetaWindow(degree_to_radian(50.0), degree_to_radian(180.0-50.0)),
    CBHemisphereGap({degree_to_radian(interval<double>::CenterWidth(0.0,40.0)),degree_to_radian(interval<double>::CenterWidth(180.0,40.0))}),
    promptrandom(ExpConfig::Setup::Get()),
    CB(Detector_t::Type_t::CB,     HistFac),
    TAPS(Detector_t::Type_t::TAPS, HistFac)
{

    h_nCaloClusters = HistFac.makeTH2D("nCaloClusters","E_{kin}^{rec} / MeV","nCaloClusters",
                                       bins_Ek, BinSettings(10),
                                       "h_nCaloClusters");
}

MCEnergyThresholds::CBTAPS_t::CBTAPS_t(Detector_t::Type_t type,
                                   const HistogramFactory& histFac) :
    Type(type),
    HistFac(Detector_t::ToString(Type), histFac, Detector_t::ToString(Type))
{

    h_Ecl_ClusterSize = HistFac.makeTH2D("E_{cl} vs. Size",  "E_{kin}^{rec} / MeV","Cluster Size",
                                         bins_Ek, bins_clSize, "h_Ecl_ClusterSize");
}

void MCEnergyThresholds::CBTAPS_t::Fill(const TCluster& caloCluster, double w) const
{
    if(caloCluster.DetectorType != Type)
        return;

    if(caloCluster.HasFlag(TCluster::Flags_t::TouchesHoleCentral))
        return;

    h_Ecl_ClusterSize->Fill(caloCluster.Energy, caloCluster.Hits.size(), w);

}

void ant::analysis::physics::MCEnergyThresholds::CBTAPS_t::Draw(canvas& c) const
{
    c << h_Ecl_ClusterSize;
}

void MCEnergyThresholds::ProcessEvent(const TEvent& event, manager_t&)
{
    triggersimu.ProcessEvent(event);

    if(!triggersimu.HasTriggered())
        return;

    const auto& data = event.Reconstructed();

    // prompt-random subtraction is easily done,
    // as clusters don't depend on tagger hits
    double taggWsum = 0.0;
    for(const TTaggerHit& taggerhit : data.TaggerHits) {
        promptrandom.SetTaggerTime(triggersimu.GetCorrectedTaggerTime(taggerhit));
        taggWsum += promptrandom.FillWeight();
    }

    for(auto& cl : data.Clusters) {
        if(CBThetaWindow.Contains(cl.Position.Theta())
           && !CBHemisphereGap.Contains(cl.Position.Phi()))
            CB.Fill(cl, taggWsum);

        TAPS.Fill(cl, taggWsum);
    }
}

namespace ant {
template<typename canvas>
canvas& operator<<(canvas&& c, const MCEnergyThresholds::CBTAPS_t& cbtaps) {
    cbtaps.Draw(c);
    return c;
}
}

void MCEnergyThresholds::ShowResult()
{
    canvas(GetName())
            << drawoption("colz")
            << padoption::enable(padoption::LogZ)
            << CB
            << endr
            << TAPS
            << endc;
}





AUTO_REGISTER_PHYSICS(MCEnergyThresholds)


