#include "CB_TimeWalk.h"

#include "expconfig/ExpConfig.h"
#include "expconfig/detectors/CB.h"

#include "base/std_ext/vector.h"
#include "base/Logger.h"
#include "base/TH_ext.h"

#include "TObjArray.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::physics;

CB_TimeWalk::CB_TimeWalk(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    noFitting(opts->Get<bool>("NoFitting", true))
{
    cb_detector = ExpConfig::Setup::GetDetector<expconfig::detector::CB>();

    BinSettings bins_energy(300); // its raw energy (not calibrated)
    BinSettings bins_channels(cb_detector->GetNChannels());

    h_timewalk =
            HistFac.makeTH3D(
                "CB TimeWalk",
                "Raw Energy",
                "Time / ns",
                "Channel",
                bins_energy,
                BinSettings(300,-50,150),
                bins_channels,
                "timewalk"
                );
}

void CB_TimeWalk::ProcessEvent(const TEvent& event, manager_t&)
{

    struct hitmapping_t {
        vector<TDetectorReadHit::Value_t> Integrals;
        vector<TDetectorReadHit::Value_t> Timings;
    };

    std::map<unsigned, hitmapping_t> hits;

    for(const TDetectorReadHit& readhit : event.Reconstructed().DetectorReadHits) {
        if(readhit.DetectorType != Detector_t::Type_t::CB)
            continue;

        auto& item = hits[readhit.Channel];

        if(readhit.ChannelType == Channel_t::Type_t::Integral) {
            std_ext::concatenate(item.Integrals, readhit.Values);
        }
        else if(readhit.ChannelType == Channel_t::Type_t::Timing) {
            std_ext::concatenate(item.Timings, readhit.Values); // passed the timing window!
        }
    }

    for(const auto& it_hit : hits) {

        const auto channel = it_hit.first;
        const hitmapping_t& item = it_hit.second;

        for(auto& integral : item.Integrals) {
            for(auto& time : item.Timings) {
                h_timewalk->Fill(integral.Uncalibrated, time.Calibrated, channel);
            }
        }
    }
}

void CB_TimeWalk::ShowResult()
{
    LOG(INFO) << "Projecting timewalk spectra...";
    // all histograms are created in subfolder now
    HistogramFactory::DirStackPush dir(HistogramFactory("Projections", HistFac));

    if(!noFitting) {
        h_timewalk_fitted =
                HistFac.makeTH2D(
                    "CB Timewalk Fitted Overview",
                    h_timewalk->GetXaxis()->GetTitle(),
                    h_timewalk->GetZaxis()->GetTitle(),
                    TH_ext::getBins(h_timewalk->GetXaxis()),
                    TH_ext::getBins(h_timewalk->GetZaxis()),
                    "h_timewalk_fitted"
                    );
    }

    canvas c_ignored(GetName()+": Ignored channels");
    c_ignored << drawoption("colz");

    for(unsigned ch=0;ch<cb_detector->GetNChannels();ch++) {
        h_timewalk->GetZaxis()->SetRange(ch+1,ch+1);
        stringstream ss_name;
        ss_name << "Ch" << ch << "_yx";
        auto proj = dynamic_cast<TH2*>(h_timewalk->Project3D(ss_name.str().c_str()));

        if(proj->GetEntries()==0) {
            delete proj;
            continue;
        }

        if(cb_detector->IsIgnored(ch))
            c_ignored << proj;

        if(noFitting || cb_detector->IsIgnored(ch))
            continue;

        LOG(INFO) << "Fitting Channel=" << ch;
        TObjArray aSlices;
        proj->FitSlicesY(nullptr, 0, -1, 0, "QNR", &aSlices);
        TH1D* means = dynamic_cast<TH1D*>(aSlices.At(1));
        for(Int_t x=0;x<means->GetNbinsX()+1;x++) {
            h_timewalk_fitted->SetBinContent(x, ch+1, means->GetBinContent(x));
        }
    }
    c_ignored << endc;
    if(!noFitting)
        canvas(GetName()) << drawoption("colz") << h_timewalk_fitted << endc;
}

AUTO_REGISTER_PHYSICS(CB_TimeWalk)