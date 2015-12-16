#include "PID_Energy.h"

#include "expconfig/ExpConfig.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::data;
using namespace ant::analysis::physics;

PID_Energy::PID_Energy(const string& name, analysis::PhysOptPtr opts) :
    Physics(name, opts)
{
    const auto detector = ExpConfig::Setup::GetDetector(Detector_t::Type_t::PID);
    const auto nChannels = detector->GetNChannels();

    const BinSettings pid_channels(nChannels);
    const BinSettings pid_rawvalues(300);


    h_pedestals = HistFac.makeTH2D(
                      "PID Pedestals",
                      "Raw ADC value",
                      "#",
                      pid_rawvalues,
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

    for(unsigned ch=0;ch<nChannels;ch++) {
        stringstream ss;
        ss << "Ch" << ch;
        h_perChannel.push_back(
                    PerChannel_t(SmartHistFactory(ss.str(), HistFac, ss.str())));
    }
}

PID_Energy::PerChannel_t::PerChannel_t(SmartHistFactory HistFac)
{
    const BinSettings cb_energy(400,0,800);
    const BinSettings pid_timing(300,-300,700);
    const BinSettings pid_rawvalues(300);
    const BinSettings pid_energy(150,0,30);

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
                    BinSettings(300,0,2000),
                    "BananaRaw"
                    );

    //    BananaTiming = HistFac.makeTH3D(
    //                       "PID Banana Timing",
    //                       "CB Energy / MeV",
    //                       "PID Energy / MeV",
    //                       "PID Timing / ns",
    //                       cb_energy,
    //                       pid_energy,
    //                       pid_timing,
    //                       "BananaTiming"
    //                       );


    TDCMultiplicity = HistFac.makeTH1D("PID TDC Multiplicity", "nHits", "#", BinSettings(10), "TDCMultiplicity");
    QDCMultiplicity = HistFac.makeTH1D("PID QDC Multiplicity", "nHits", "#", BinSettings(10), "QDCMultiplicity");
}

void PID_Energy::ProcessEvent(const Event& event)
{
    // pedestals, best determined from clusters with energy information only
    for(const Cluster& cluster : event.Reconstructed.AllClusters) {
        if(!(cluster.Detector & Detector_t::Type_t::PID))
            continue;

        // only consider one cluster PID hits
        if(cluster.Hits.size() != 1)
            continue;
        const Cluster::Hit& pidhit = cluster.Hits.front();

        PerChannel_t& h = h_perChannel[pidhit.Channel];

        double pedestal = numeric_limits<double>::quiet_NaN();
        double timing = numeric_limits<double>::quiet_NaN();

        unsigned nPedestals = 0;
        unsigned nTimings = 0;
        for(const Cluster::Hit::Datum& datum : pidhit.Data) {
            if(datum.Type == Channel_t::Type_t::Pedestal) {
                pedestal = datum.Value;
                nPedestals++;
            }
            if(datum.Type == Channel_t::Type_t::Timing) {
                timing = datum.Value;
                nTimings++;
            }
        }

        h.QDCMultiplicity->Fill(nPedestals);
        h.TDCMultiplicity->Fill(nTimings);

        if(!std::isfinite(pedestal))
            continue;
        if(nPedestals > 1)
            continue;
        if(nTimings > 1)
            continue;

        h_pedestals->Fill(pedestal, pidhit.Channel);

        if(std::isfinite(timing))
            h.PedestalTiming->Fill(timing, pedestal);
        else
            h.PedestalNoTiming->Fill(pedestal);

    }


    // bananas per channel histograms
    for(const auto& candidate : event.Reconstructed.Candidates) {
        // only candidates with one cluster in CB and one cluster in PID
        if(candidate->Clusters.size() != 2)
            continue;
        const bool cb_and_pid = candidate->Detector & Detector_t::Type_t::CB &&
                                candidate->Detector & Detector_t::Type_t::PID;
        if(!cb_and_pid)
            continue;

        // search for PID cluster
        auto pid_cluster = candidate->FindFirstCluster(Detector_t::Type_t::PID);

        h_bananas->Fill(candidate->ClusterEnergy,
                        candidate->VetoEnergy,
                        pid_cluster->CentralElement);

        // per channel histograms
        PerChannel_t& h = h_perChannel[pid_cluster->CentralElement];

        // fill the banana
        h.Banana->Fill(candidate->ClusterEnergy,
                       candidate->VetoEnergy);


        double pedestal = numeric_limits<double>::quiet_NaN();
        //double timing = numeric_limits<double>::quiet_NaN();
        for(const Cluster::Hit& pidhit : pid_cluster->Hits) {
            for(const Cluster::Hit::Datum& datum : pidhit.Data) {
                if(datum.Type == Channel_t::Type_t::Pedestal) {
                    pedestal = datum.Value;
                }
                //                if(datum.Type == Channel_t::Type_t::Timing) {
                //                    timing = datum.Value;
                //                }
            }
        }

        h.BananaRaw->Fill(candidate->ClusterEnergy, pedestal);
        //h.BananaTiming->Fill(candidate->ClusterEnergy(), candidate->VetoEnergy(), timing);

    }
}

void PID_Energy::ShowResult()
{
    canvas(GetName())
            << drawoption("colz") << h_pedestals
            << endc;
}

AUTO_REGISTER_PHYSICS(PID_Energy)

