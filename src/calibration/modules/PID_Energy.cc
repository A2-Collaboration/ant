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


    h_pedestals = HistFac.makeTH2D(
                      "PID Pedestals",
                      "Raw ADC value",
                      "#",
                      pid_rawvalues,
                      pid_channels,
                      "Pedestals");

    for(unsigned ch=0;ch<nChannels;ch++) {
        stringstream ss;
        ss << "Ch" << ch;
        h_perChannel.push_back(
                    PerChannel_t(SmartHistFactory(ss.str(), HistFac, ss.str())));
    }

}

PID_Energy::ThePhysics::PerChannel_t::PerChannel_t(SmartHistFactory HistFac)
{
    const BinSettings cb_energy(200,0,800);
    const BinSettings pid_timing(100,-100,100);
    const BinSettings pid_rawvalues(300);
    const BinSettings pid_energy(50,0,30);

    PedestalTiming = HistFac.makeTH2D(
                         "PID Pedestals Timing",
                         "Timing / ns",
                         "Raw ADC value",
                         pid_timing,
                         pid_rawvalues,
                         "PedestalTiming");

    PedestalNoTiming = HistFac.makeTH1D(
                           "PID Pedestals No Timing",
                           "Raw ADC value",
                           "#",
                           pid_rawvalues,
                           "PedestalNoTiming");

    Banana = HistFac.makeTH2D(
                 "PID Banana",
                 "CB Energy / MeV",
                 "PID Energy / MeV",
                 cb_energy,
                 pid_energy,
                 "Banana"
                );

    BananaRaw = HistFac.makeTH2D(
                "PID Banana Raw",
                "CB Energy / MeV",
                "PID ADC Value",
                cb_energy,
                pid_rawvalues,
                "BananaRaw"
                );

    BananaTiming = HistFac.makeTH3D(
                       "PID Banana Timing",
                       "CB Energy / MeV",
                       "PID Energy / MeV",
                       "PID Timing / ns",
                       cb_energy,
                       pid_energy,
                       pid_timing,
                       "BananaTiming"
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
        const Cluster::Hit& pidhit = cluster.Hits.front();

        PerChannel_t& h = h_perChannel[pidhit.Channel];

        double pedestal = numeric_limits<double>::quiet_NaN();
        double timing = numeric_limits<double>::quiet_NaN();

        for(const Cluster::Hit::Datum& datum : pidhit.Data) {
            if(datum.Type == Channel_t::Type_t::Pedestal) {
                pedestal = datum.Value;
            }
            if(datum.Type == Channel_t::Type_t::Timing) {
                timing = datum.Value;
            }
        }

        if(!std::isfinite(pedestal))
            continue;


        h_pedestals->Fill(pedestal, pidhit.Channel);

        if(std::isfinite(timing))
            h.PedestalTiming->Fill(timing, pedestal);
        else
            h.PedestalNoTiming->Fill(pedestal);

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

        PerChannel_t& h = h_perChannel[pid_cluster->CentralElement];

        // fill the banana
        h.Banana->Fill(candidate->ClusterEnergy(),
                       candidate->VetoEnergy());


        double pedestal = numeric_limits<double>::quiet_NaN();
        double timing = numeric_limits<double>::quiet_NaN();
        for(const Cluster::Hit& pidhit : pid_cluster->Hits) {
            for(const Cluster::Hit::Datum& datum : pidhit.Data) {
                if(datum.Type == Channel_t::Type_t::Pedestal) {
                    pedestal = datum.Value;
                }
                if(datum.Type == Channel_t::Type_t::Timing) {
                    timing = datum.Value;
                }
            }
        }

        h.BananaRaw->Fill(candidate->ClusterEnergy(), pedestal);
        h.BananaTiming->Fill(candidate->ClusterEnergy(), candidate->VetoEnergy(), timing);

    }
}

void PID_Energy::ThePhysics::Finish()
{
}

void PID_Energy::ThePhysics::ShowResult()
{
    canvas(GetName())
            << drawoption("colz") << h_pedestals
            << endc;
}





