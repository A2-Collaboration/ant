#include "PID_Energy.h"

#include "analysis/plot/HistogramFactories.h"
#include "analysis/data/Event.h"

#include "expconfig/detectors/PID.h"

#include "tree/TDataRecord.h"

#include "base/Logger.h"

#include <list>


using namespace std;
using namespace ant;
using namespace ant::calibration;
using namespace ant::analysis;
using namespace ant::analysis::data;

PID_Energy::PID_Energy(
        std::shared_ptr<expconfig::detector::PID> pid,
        std::shared_ptr<DataManager> calmgr,
        Calibration::Converter::ptr_t converter,
        double defaultPedestal,
        double defaultGain,
        double defaultThreshold,
        double defaultRelativeGain
        ) :
    Energy(pid->Type,
           calmgr,
           converter,
           defaultPedestal,
           defaultGain,
           defaultThreshold,
           defaultRelativeGain),
    pid_detector(pid)
{

}


PID_Energy::~PID_Energy()
{

}

unique_ptr<analysis::Physics> PID_Energy::GetPhysicsModule()
{
    return std_ext::make_unique<ThePhysics>(GetName(),
                                            pid_detector->GetNChannels());
}

void PID_Energy::GetGUIs(std::list<std::unique_ptr<gui::Manager_traits> >& guis)
{
    guis.emplace_back(std_ext::make_unique<GUI_Pedestals>(
                          GetName(),
                          Pedestals,
                          calibrationManager,
                          pid_detector
                          ));
}

PID_Energy::ThePhysics::ThePhysics(const string& name,
                                   unsigned nChannels):
    Physics(name)
{
    const BinSettings pid_channels(nChannels);

    h_pedestals = HistFac.makeTH2D(
                      "PID Pedestals",
                      "Raw ADC value",
                      "#",
                      BinSettings(300),
                      pid_channels,
                      "Pedestals");
    h_bananas =
            HistFac.makeTH3D(
                "PID Bananas",
                "CB Energy / MeV",
                "PID Energy / MeV",
                "Channel",
                BinSettings(400,0,800),
                BinSettings(100,0,30),
                pid_channels,
                "Bananas"
                );
}

void PID_Energy::ThePhysics::ProcessEvent(const data::Event& event)
{

    // pedestals, best determined from clusters with energy information only
    for(const Cluster& cluster : event.Reconstructed().AllClusters()) {
        if(!(cluster.Detector & Detector_t::Type_t::PID))
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
    for(const auto& candidate : event.Reconstructed().Candidates()) {
        if(candidate->Clusters.size() != 2)
            continue;
        if(candidate->Detector() & Detector_t::Type_t::CB &&
           candidate->Detector() & Detector_t::Type_t::PID
           )
        {
            // search for PID cluster
            auto pid_cluster = candidate->FindFirstCluster(Detector_t::Type_t::PID);
            h_bananas->Fill(candidate->ClusterEnergy(),
                            candidate->VetoEnergy(),
                            pid_cluster->CentralElement);
        }
    }
}

void PID_Energy::ThePhysics::Finish()
{
}

void PID_Energy::ThePhysics::ShowResult()
{
    canvas(GetName())
            << drawoption("colz") << h_pedestals
            << drawoption("colz") << h_bananas->Project3D("zy")
            << endc;
}


