#include "TAPS_ShortEnergy.h"

#include "expconfig/ExpConfig.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::physics;



TAPS_ShortEnergy::TAPS_ShortEnergy(const string& name, OptionsPtr opts) :
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
                      BinSettings(400,0,4),
                      taps_channels,
                      "rel_gamma");
}

void TAPS_ShortEnergy::ProcessEvent(const TEvent& event, manager_t&)
{
    // pedestals
    for(const TDetectorReadHit& readhit : event.Reconstructed().DetectorReadHits) {
        if(readhit.DetectorType != Detector_t::Type_t::TAPS)
            continue;
        if(readhit.ChannelType != Channel_t::Type_t::IntegralShort)
            continue;
        /// \todo check for timing hit?
        /// \todo check for trigger pattern?
        for(const double& value : readhit.Converted)
            h_pedestals->Fill(value, readhit.Channel);
    }

    for(const TCandidatePtr& c : event.Reconstructed().Candidates) {
        if(c->Detector & Detector_t::Any_t::TAPS_Apparatus
           && c->VetoEnergy < 0.5)
        {
            const auto& cluster = c->FindCaloCluster();

            if(cluster)
                for(const TClusterHit& clusterhit : cluster->Hits) {

                    if(clusterhit.Channel == cluster->CentralElement) {
                        double central_e = numeric_limits<double>::quiet_NaN();
                        double short_e = numeric_limits<double>::quiet_NaN();
                        for(const TClusterHit::Datum& datum : clusterhit.Data) {

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