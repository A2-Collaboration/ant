#include "TAPSVeto_Energy.h"
#include "analysis/plot/HistogramFactories.h"
#include "analysis/data/Event.h"
#include "analysis/utils/combinatorics.h"

#include "tree/TDataRecord.h"

#include "expconfig/detectors/TAPSVeto.h"

#include <list>


using namespace std;
using namespace ant;
using namespace ant::calibration;
using namespace ant::analysis;
using namespace ant::analysis::data;

TAPSVeto_Energy::TAPSVeto_Energy(std::shared_ptr<expconfig::detector::TAPSVeto> tapsveto,
                                 std::shared_ptr<DataManager> calmgr,
                                 Calibration::Converter::ptr_t converter,
                                 double defaultPedestal,
                                 double defaultGain,
                                 double defaultThreshold,
                                 double defaultRelativeGain):
    Energy(Detector_t::Type_t::TAPSVeto,
           calmgr,
           converter,
           defaultPedestal,
           defaultGain,
           defaultThreshold,
           defaultRelativeGain),
    tapsveto_detector(tapsveto)
{

}

TAPSVeto_Energy::ThePhysics::ThePhysics(const string& name, unsigned nChannels):
    Physics(name)
{
    const BinSettings tapsveto_channels(nChannels);

    h_pedestals = HistFac.makeTH2D(
                      "TAPSVeto Pedestals",
                      "Raw ADC value",
                      "#",
                      BinSettings(300),
                      tapsveto_channels,
                      "Pedestals");
    h_bananas =
            HistFac.makeTH3D(
                "TAPSVeto Bananas",
                "TAPS LG Energy / MeV",
                "TAPSVeto Energy / MeV",
                "Channel",
                BinSettings(400,0,800),
                BinSettings(100,0,30),
                tapsveto_channels,
                "Bananas"
                );}

void TAPSVeto_Energy::ThePhysics::ProcessEvent(const Event& event)
{
    const auto& cands = event.Reconstructed().Candidates();

    // pedestals
    for(const Cluster& cluster : event.Reconstructed().AllClusters()) {
        if(!(cluster.Detector & Detector_t::Type_t::TAPSVeto))
            continue;
        for(const Cluster::Hit& clusterhit : cluster.Hits) {
            /// \todo check for timing hit?
            for(const Cluster::Hit::Datum& datum : clusterhit.Data) {
                if(datum.Type != Channel_t::Type_t::Pedestal)
                    continue;
                h_pedestals->Fill(datum.Value, clusterhit.Channel);
            }
        }
    }

    // bananas
    for(const auto& candidate : cands) {
        if(candidate->Clusters.size() != 2)
            continue;
        if(candidate->Detector() & Detector_t::Type_t::TAPS &&
           candidate->Detector() & Detector_t::Type_t::TAPSVeto
           )
        {
            // search for TAPSVeto cluster
            auto tapsveto_cluster = candidate->FindFirstCluster(Detector_t::Type_t::TAPSVeto);
            h_bananas->Fill(candidate->ClusterEnergy(),
                            candidate->VetoEnergy(),
                            tapsveto_cluster->CentralElement);
        }
    }

}

void TAPSVeto_Energy::ThePhysics::Finish()
{
}

void TAPSVeto_Energy::ThePhysics::ShowResult()
{
    canvas(GetName())
            << drawoption("colz") << h_pedestals
            << drawoption("colz") << h_bananas->Project3D("zy")
            << endc;
}

unique_ptr<analysis::Physics> TAPSVeto_Energy::GetPhysicsModule()
{
    return std_ext::make_unique<ThePhysics>(GetName(), tapsveto_detector->GetNChannels());
}

void TAPSVeto_Energy::GetGUIs(std::list<std::unique_ptr<gui::Manager_traits> >& guis) {
    guis.emplace_back(std_ext::make_unique<GUI_Pedestals>(
                          GetName(),
                          Pedestals,
                          calibrationManager,
                          tapsveto_detector
                          ));
}
