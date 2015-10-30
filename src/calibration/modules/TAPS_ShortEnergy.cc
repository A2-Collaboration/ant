#include "TAPS_ShortEnergy.h"

#include "expconfig/detectors/TAPS.h"

#include "analysis/plot/HistogramFactories.h"
#include "analysis/data/Event.h"
#include "analysis/utils/combinatorics.h"

#include "tree/TDataRecord.h"

#include <list>
#include <cmath>

using namespace std;
using namespace ant;
using namespace ant::calibration;
using namespace ant::analysis;
using namespace ant::analysis::data;

TAPS_ShortEnergy::TAPS_ShortEnergy(std::shared_ptr<expconfig::detector::TAPS> taps,
        std::shared_ptr<DataManager> calmgr,
        Calibration::Converter::ptr_t converter,
        double defaultPedestal,
        double defaultGain,
        double defaultThreshold,
        double defaultRelativeGain) :
    Energy(Detector_t::Type_t::TAPS,
           calmgr,
           converter,
           defaultPedestal,
           defaultGain,
           defaultThreshold,
           defaultRelativeGain,
           Channel_t::Type_t::IntegralShort),
    taps_detector(taps)
{

}

TAPS_ShortEnergy::ThePhysics::ThePhysics(const string& name, shared_ptr<expconfig::detector::TAPS> taps) :
    Physics(name),
    taps_detector(taps)
{
    const BinSettings taps_channels(taps->GetNChannels());

    h_pedestals = HistFac.makeTH2D(
                      "TAPS ShortGate Pedestals",
                      "Raw ADC value",
                      "#",
                      BinSettings(300),
                      taps_channels,
                      "Pedestals");
}

void TAPS_ShortEnergy::ThePhysics::ProcessEvent(const Event& event)
{
    // pedestals
    for(const Cluster& cluster : event.Reconstructed().AllClusters()) {
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
}

void TAPS_ShortEnergy::ThePhysics::Finish()
{
}

void TAPS_ShortEnergy::ThePhysics::ShowResult()
{
    canvas(GetName()) << drawoption("colz") << h_pedestals
                      << endc;
}

unique_ptr<analysis::Physics> TAPS_ShortEnergy::GetPhysicsModule()
{
    return std_ext::make_unique<ThePhysics>(GetName(), taps_detector);
}

void TAPS_ShortEnergy::GetGUIs(std::list<std::unique_ptr<gui::Manager_traits> >& guis)
{
    guis.emplace_back(std_ext::make_unique<GUI_Pedestals>(
                          GetName(),
                          Pedestals,
                          calibrationManager,
                          taps_detector
                          ));
}
