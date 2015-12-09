#include "CB_TimeWalk.h"

#include "expconfig/ExpConfig.h"
#include "expconfig/detectors/CB.h"

#include "base/Logger.h"

#include "TObjArray.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::data;
using namespace ant::analysis::physics;

CB_TimeWalk::CB_TimeWalk(const string& name, analysis::PhysOptPtr opts) :
    Physics(name, opts)
{
    cb_detector = ExpConfig::Setup::GetDetector<expconfig::detector::CB>();

    BinSettings bins_energy(400,0,500);
    BinSettings bins_channels(cb_detector->GetNChannels());

    h_timewalk =
            HistFac.makeTH3D(
                "CB TimeWalk",
                "Energy / MeV",
                "Time / ns",
                "Channel",
                bins_energy,
                BinSettings(100,-100,100),
                bins_channels,
                "timewalk"
                );
    h_timewalk_overview =
            HistFac.makeTH2D(
                "CB Timewalk Overview",
                "Energy / MeV",
                "Channel",
                bins_energy,
                bins_channels,
                "timewalk_overview"
                );
}

void CB_TimeWalk::ProcessEvent(const analysis::data::Event& event)
{
    for(const auto& cand: event.Reconstructed.Candidates) {
        for(const Cluster& cluster: cand->Clusters) {
            if(cluster.Detector != Detector_t::Type_t::CB)
                continue;
            for(const Cluster::Hit& hit : cluster.Hits) {
                // found the hit of the central element
                // now search for its timing information
                double time = numeric_limits<double>::quiet_NaN();
                double energy = numeric_limits<double>::quiet_NaN();
                for(const Cluster::Hit::Datum& d : hit.Data) {
                    if(d.Type == Channel_t::Type_t::Timing)
                        time = d.Value;
                    if(d.Type == Channel_t::Type_t::Integral)
                        energy = d.Value;
                }
                h_timewalk->Fill(energy, time, hit.Channel);
            }
        }
    }
}

void CB_TimeWalk::ShowResult()
{
    for(unsigned ch=0;ch<cb_detector->GetNChannels();ch++) {
        if(cb_detector->IsIgnored(ch))
            continue;
        LOG(INFO) << "Fitting Channel=" << ch;
        h_timewalk->GetZaxis()->SetRange(ch,ch+1);
        stringstream ss_name;
        ss_name << "Ch" << ch << "_yx";
        TH2* proj = dynamic_cast<TH2*>(h_timewalk->Project3D(ss_name.str().c_str()));
        TObjArray aSlices;
        proj->FitSlicesY(nullptr, 0, -1, 0, "QNR", &aSlices);
        TH1D* means = dynamic_cast<TH1D*>(aSlices.At(1));
        for(Int_t x=0;x<means->GetNbinsX()+1;x++) {
            h_timewalk_overview->SetBinContent(x, ch+1, means->GetBinContent(x));
        }
        delete proj;
    }
    canvas(GetName()) << drawoption("colz") << h_timewalk_overview << endc;
}

AUTO_REGISTER_PHYSICS(CB_TimeWalk)