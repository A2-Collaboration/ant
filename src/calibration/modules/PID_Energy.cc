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
    const BinSettings pid_rawvalues(300);
    const BinSettings cb_energy(400,0,800);

    h_pedestals = HistFac.makeTH2D(
                      "PID Pedestals",
                      "Raw ADC value",
                      "#",
                      pid_rawvalues,
                      pid_channels,
                      "Pedestals");
    h_pedestals_timing = HistFac.makeTH3D(
                             "PID Pedestals",
                             "Raw ADC value",
                             "Timing / ns",
                             "#",
                             pid_rawvalues,
                             BinSettings(100, -50, 50),
                             pid_channels,
                             "Pedestals_timing");

    h_bananas =
            HistFac.makeTH3D(
                "PID Bananas",
                "CB Energy / MeV",
                "PID Energy / MeV",
                "Channel",
                cb_energy,
                BinSettings(100,0,30),
                pid_channels,
                "Bananas"
                );
    h_bananas_raw =
            HistFac.makeTH3D(
                "PID Bananas Raw",
                "CB Energy / MeV",
                "PID ADC Value",
                "Channel",
                BinSettings(400,0,800),
                pid_rawvalues,
                pid_channels,
                "Bananas_raw"
                );
}

void PID_Energy::ThePhysics::ProcessEvent(const data::Event& event)
{
    // pedestals, best determined from clusters with energy information only
    for(const Cluster& cluster : event.Reconstructed().AllClusters()) {
        if(!(cluster.Detector & Detector_t::Type_t::PID))
            continue;

        // only consider one cluster PID hits
        if(cluster.Hits.size() != 1)
            continue;
        const Cluster::Hit& clusterhit = cluster.Hits.front();

        double pedestal = numeric_limits<double>::quiet_NaN();
        double timing = numeric_limits<double>::quiet_NaN();

        for(const Cluster::Hit::Datum& datum : clusterhit.Data) {
            if(datum.Type == Channel_t::Type_t::Pedestal) {
                pedestal = datum.Value;
            }
            if(datum.Type == Channel_t::Type_t::Timing) {
                timing = datum.Value;
            }
        }

        h_pedestals->Fill(pedestal, clusterhit.Channel);
        h_pedestals_timing->Fill(pedestal, timing, clusterhit.Channel);

    }

    // bananas
    for(const auto& candidate : event.Reconstructed().Candidates()) {
        // only candidates with one cluster in CB and one cluster in PID
        if(candidate->Clusters.size() != 2)
            continue;
        const bool cb_and_pid = candidate->Detector() & Detector_t::Type_t::CB &&
                                candidate->Detector() & Detector_t::Type_t::PID;
        if(!cb_and_pid)
            continue;

        // search for PID cluster
        auto pid_cluster = candidate->FindFirstCluster(Detector_t::Type_t::PID);

        // fill the banana
        h_bananas->Fill(candidate->ClusterEnergy(),
                        candidate->VetoEnergy(),
                        pid_cluster->CentralElement);


        for(const Cluster::Hit& pidhit : pid_cluster->Hits) {
            for(const Cluster::Hit::Datum& datum : pidhit.Data) {
                if(datum.Type != Channel_t::Type_t::Pedestal)
                    continue;
                h_bananas_raw->Fill(candidate->ClusterEnergy(),
                                    datum.Value,
                                    pid_cluster->CentralElement);
            }
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
            << drawoption("colz") << h_pedestals_timing->Project3D("zy")
            << drawoption("colz") << h_bananas->Project3D("zy")
            << drawoption("colz") << h_bananas_raw->Project3D("zy")
            << endc;
}


