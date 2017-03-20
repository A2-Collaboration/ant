#include "CB_TimeWalk.h"

#include "expconfig/ExpConfig.h"
#include "expconfig/detectors/CB.h"
#include "calibration/converters/CATCH_TDC.h"

#include "base/std_ext/container.h"
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

    // logarithmic from about 10 (raw units, threshold is at 16)
    // to 2^14, which is maximum ADC value
    // we need to full energy range for good TimeWalk correction!
    const BinSettings bins_energy(500,std::log10(10),std::log10(1 << 14));

    const BinSettings bins_channels(cb_detector->GetNChannels());
    const auto bins_timings = BinSettings::RoundToBinSize({200,-50,150}, calibration::converter::Gains::CATCH_TDC);
    h_timewalk =
            HistFac.makeTH3D(
                "CB TimeWalk",
                "log_{10}(RawEnergy)",
                "Time / ns",
                "Channel",
                bins_energy,
                bins_timings,
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
                h_timewalk->Fill(std::log10(integral.Uncalibrated), time.Calibrated, channel);
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

    canvas c(GetName());

    h_timewalk->GetZaxis()->SetRange(1,cb_detector->GetNChannels());
    c << drawoption("colz") << dynamic_cast<TH2*>(h_timewalk->Project3D("yx"));

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
    if(!noFitting)
        c << h_timewalk_fitted;
    c << endc;
    c_ignored << endc;

}

AUTO_REGISTER_PHYSICS(CB_TimeWalk)