#include "PID_Energy.h"

#include "expconfig/ExpConfig.h"

#include "base/std_ext/vector.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::physics;

PID_Energy::PID_Energy(const string& name, OptionsPtr opts) :
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
                BinSettings(200,0,18),
                pid_channels,
                "Bananas"
                );

    for(unsigned ch=0;ch<nChannels;ch++) {
        stringstream ss;
        ss << "Ch" << ch;
        h_perChannel.push_back(
                    PerChannel_t(HistogramFactory(ss.str(), HistFac, ss.str())));
    }
}

PID_Energy::PerChannel_t::PerChannel_t(HistogramFactory HistFac)
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

void PID_Energy::ProcessEvent(const TEvent& event, manager_t&)
{
    // pedestals, best determined from clusters with energy information only

    struct hitmapping_t {
        vector<double> Pedestals;
        vector<double> Timings;
    };

    std::map<unsigned, hitmapping_t> hits;

    for(const TDetectorReadHit& readhit : event.Reconstructed().DetectorReadHits) {
        if(readhit.DetectorType != Detector_t::Type_t::PID)
            continue;

        auto& item = hits[readhit.Channel];

        if(readhit.ChannelType == Channel_t::Type_t::Integral) {
            std_ext::concatenate(item.Pedestals, readhit.Converted);
        }
        else if(readhit.ChannelType == Channel_t::Type_t::Timing) {
            std_ext::concatenate(item.Timings, readhit.Values); // passed the timing window!
        }
    }

    for(const auto& it_hit : hits) {

        const auto channel = it_hit.first;
        const hitmapping_t& item = it_hit.second;

        PerChannel_t& h = h_perChannel[channel];

        h.QDCMultiplicity->Fill(item.Pedestals.size());
        h.TDCMultiplicity->Fill(item.Timings.size());


        if(item.Pedestals.size() != 1)
            continue;
        if(item.Timings.size()>1)
            continue;

        const auto& pedestal = item.Pedestals.front();

        h_pedestals->Fill(pedestal, channel);

        if(item.Timings.size()==1)
            h.PedestalTiming->Fill(item.Timings.front(), pedestal);
        else
            h.PedestalNoTiming->Fill(pedestal);

    }

    // bananas per channel histograms
    for(const auto& candidate : event.Reconstructed().Candidates) {
        // only candidates with one cluster in CB and one cluster in PID
        if(candidate->Clusters.size() != 2)
            continue;
        const bool cb_and_pid = candidate->Detector & Detector_t::Type_t::CB &&
                                candidate->Detector & Detector_t::Type_t::PID;
        if(!cb_and_pid)
            continue;

        // search for PID cluster
        const TClusterPtr& pid_cluster = candidate->FindFirstCluster(Detector_t::Type_t::PID);

        h_bananas->Fill(candidate->CaloEnergy,
                        candidate->VetoEnergy,
                        pid_cluster->CentralElement);

        // per channel histograms
        PerChannel_t& h = h_perChannel[pid_cluster->CentralElement];

        // fill the banana
        h.Banana->Fill(candidate->CaloEnergy,
                       candidate->VetoEnergy);

        // is there an pedestal available?
        const auto it_hit = hits.find(pid_cluster->CentralElement);
        if(it_hit == hits.end()) {
            continue;
        }
        const auto& pedestals = it_hit->second.Pedestals;
        if(pedestals.size() != 1)
            continue;

        const auto& pedestal = pedestals.front();

        h.BananaRaw->Fill(candidate->CaloEnergy, pedestal);
        //h.BananaTiming->Fill(candidate->ClusterEnergy(), candidate->VetoEnergy(), timing);

    }
}

void PID_Energy::ShowResult()
{
    canvas(GetName())
            << drawoption("colz") << h_pedestals
            << endc;
    canvas c_bananas(GetName()+": Bananas");
    c_bananas << drawoption("colz");
    for(auto& h : h_perChannel)
        c_bananas << h.Banana;
    c_bananas << endc;
}

AUTO_REGISTER_PHYSICS(PID_Energy)

