#include "TAPS_ShortEnergy.h"

#include "expconfig/ExpConfig.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::data;
using namespace ant::analysis::physics;



TAPS_ShortEnergy::TAPS_ShortEnergy(const string& name, analysis::PhysOptPtr opts) :
    Physics(name, opts)
{
    auto taps = ExpConfig::Setup::GetDetector(Detector_t::Type_t::TAPS);

    const BinSettings taps_channels(taps->GetNChannels());

    h_pedestals = HistFac.makeTH2D(
                      "TAPS ShortGate Pedestals",
                      "Raw ADC value",
                      "#",
                      BinSettings(300),
                      taps_channels,
                      "Pedestals");

    h_rel_gamma = HistFac.makeTH2D(
                      "TAPS E_{S} / E_{L}",
                      "rel",
                      "#",
                      BinSettings(400,-10,10),
                      taps_channels,
                      "RelativeGains");
}

void TAPS_ShortEnergy::ProcessEvent(const Event& event)
{
    // pedestals
    for(const Cluster& cluster : event.Reconstructed.AllClusters) {
        if(!(cluster.Detector & Detector_t::Type_t::TAPS))
            continue;
        for(const Cluster::Hit& clusterhit : cluster.Hits) {
            /// \todo check for timing hit?
            /// \todo check for trigger pattern?
            for(const Cluster::Hit::Datum& datum : clusterhit.Data) {

                if(datum.Type != Channel_t::Type_t::PedestalShort)
                    continue;
                h_pedestals->Fill(datum.Value, clusterhit.Channel);

            }
        }
    }

    for(const auto& c : event.Reconstructed.Candidates) {
        if(c->VetoEnergy < 0.5) {
            const auto& cluster = c->FindCaloCluster();

            if(cluster)
                for(const Cluster::Hit& clusterhit : cluster->Hits) {

                    if(clusterhit.Channel == cluster->CentralElement) {
                        double central_e = numeric_limits<double>::quiet_NaN();
                        double short_e = numeric_limits<double>::quiet_NaN();
                        for(const Cluster::Hit::Datum& datum : clusterhit.Data) {

                            if(datum.Type == Channel_t::Type_t::Integral)
                                central_e = datum.Value;

                            if(datum.Type == Channel_t::Type_t::IntegralShort)
                                short_e = datum.Value;
                        }

                        h_rel_gamma->Fill(short_e / central_e, clusterhit.Channel);

                    }
                }
        }
    }
}

void TAPS_ShortEnergy::ShowResult()
{
    canvas(GetName()) << drawoption("colz") << h_pedestals
                      << h_rel_gamma
                      << endc;
}

AUTO_REGISTER_PHYSICS(TAPS_ShortEnergy)