#include "Thresholds.h"

#include "expconfig/ExpConfig.h"
#include "root-addons/cbtaps_display/TH2CB.h"
#include "base/std_ext/string.h"
#include "base/TH_ext.h"

using namespace std;
using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::physics;

Thresholds::Thresholds(const Detector_t::Type_t& detectorType,
                       const BinSettings& bins_x,
                       const string& name, OptionsPtr opts) :
    Physics(name, opts),
    Detector(ExpConfig::Setup::GetDetector(detectorType)),
    badTDCthreshold(opts->Get<double>("BadTDCThreshold", 7.0))
{
    const BinSettings bins_ch(Detector->GetNChannels());

    hThresholds_Raw = HistFac.makeTH2D("Thresholds Raw","Raw","Element",
                                       BinSettings(300), bins_ch,
                                       "hThresholds_Raw");

    hThresholds_ADC = HistFac.makeTH2D("Thresholds ADC","Energy","Element",
                                       bins_x, bins_ch,
                                       "hThresholds_ADC");

    hThresholds_TDC = HistFac.makeTH2D("Thresholds TDC","Energy","Element",
                                       bins_x, bins_ch,
                                       "hThresholds_TDC");

    hThresholds_Raw_TDC= HistFac.makeTH2D("Thresholds Raw TDC","Raw","Element",
                                          BinSettings(300), bins_ch,
                                          "hThresholds_Raw_TDC");

    // makes only sense if Ant is running with '-S IncludeIgnoredElements=1'
    hADCnoTDC = HistFac.makeTH1D(
                formatter() << "ADC>" << badTDCthreshold << "MeV but no TDC","Channel","",
                                      bins_ch, "hADCnoTDC");

    hADC      = HistFac.makeTH1D(
                formatter() << "ADC>" << badTDCthreshold << "MeV", "Channel","",
                                      bins_ch, "hADC");
    hADCnoTDC_norm  = HistFac.makeTH1D(
                formatter() << "ADC>" << badTDCthreshold << "MeV Normalized", "Channel","",
                                      bins_ch, "hADCnoTDC_norm");
}

void Thresholds::ProcessEvent(const TEvent& event, manager_t&)
{
    struct hit_t {
        double Energy = std_ext::NaN;
        double Time = std_ext::NaN;
        double Pedestal = std_ext::NaN;
    };

    map<unsigned, hit_t> hits;
    for(const TDetectorReadHit& readhit : event.Reconstructed().DetectorReadHits) {
        if(readhit.DetectorType != Detector->Type)
            continue;
        if(readhit.Values.empty())
            continue;
        auto& hit = hits[readhit.Channel];
        if(readhit.ChannelType == Channel_t::Type_t::Integral) {
            hit.Pedestal = readhit.Values.front().Uncalibrated;
            hit.Energy = readhit.Values.front().Calibrated;
        }
        if(readhit.ChannelType == Channel_t::Type_t::Timing)
            hit.Time = readhit.Values.front().Calibrated;
    }
    for(const auto& it_hit : hits) {
        auto ch = it_hit.first;
        const hit_t& hit = it_hit.second;
        hThresholds_Raw->Fill(hit.Pedestal, ch);
        hThresholds_ADC->Fill(hit.Energy, ch);

        if(hit.Energy > badTDCthreshold)
            hADC->Fill(ch);

        if(isfinite(hit.Time)){
            hThresholds_TDC->Fill(hit.Energy, ch);
            hThresholds_Raw_TDC->Fill(hit.Pedestal,ch);
        }
        else if(hit.Energy > badTDCthreshold)
            hADCnoTDC->Fill(ch);
    }
}

void Thresholds::ShowResult()
{
    auto CBPlot = [&] (TH1* h) {
        auto cbhist = new TH2CB(Form("cbhist_%s",h->GetName()), h->GetTitle());
        cbhist->FillElements(*h);
        for(unsigned ch=0;ch<Detector->GetNChannels();ch++) {
            if(Detector->IsIgnored(ch)) {
                cbhist->CreateMarker(ch);
            }
        }
        return cbhist;
    };

    canvas c(GetName());
    c << drawoption("colz")
      << hThresholds_Raw
      << hThresholds_ADC
      << hThresholds_TDC
      << hADCnoTDC;
    if(Detector->Type == Detector_t::Type_t::CB) {
        c << CBPlot(hADCnoTDC);
        c << CBPlot(hADCnoTDC_norm);
    }
    c << hADCnoTDC_norm;
    c << endc;
}

void Thresholds::Finish()
{
    for(int i=1;i<hADCnoTDC->GetNbinsX();++i) {
        const auto a = hADCnoTDC->GetBinContent(i);
        const auto b = hADC->GetBinContent(i);
        hADCnoTDC_norm->SetBinContent(i, b!=0.0 ? a/b : 0.0 );
    }
}

struct EPT_Thresholds : Thresholds {
    EPT_Thresholds(const std::string& name, OptionsPtr opts) :
        Thresholds(Detector_t::Type_t::EPT,
                   BinSettings(2500),
                   name, opts)
    {}
};
AUTO_REGISTER_PHYSICS(EPT_Thresholds)

struct CB_Thresholds : Thresholds {
    CB_Thresholds(const std::string& name, OptionsPtr opts) :
        Thresholds(Detector_t::Type_t::CB,
                   BinSettings(300,0,20),
                   name, opts)
    {}
};
AUTO_REGISTER_PHYSICS(CB_Thresholds)

struct PID_Thresholds : Thresholds {
    PID_Thresholds(const std::string& name, OptionsPtr opts) :
        Thresholds(Detector_t::Type_t::PID,
                   BinSettings(300,0,20),
                   name, opts)
    {}
};
AUTO_REGISTER_PHYSICS(PID_Thresholds)

struct TAPS_Thresholds : Thresholds {
    TAPS_Thresholds(const std::string& name, OptionsPtr opts) :
        Thresholds(Detector_t::Type_t::TAPS,
                   BinSettings(100,0,20),
                   name, opts)
    {}
};
AUTO_REGISTER_PHYSICS(TAPS_Thresholds)

struct TAPSVeto_Thresholds : Thresholds {
    TAPSVeto_Thresholds(const std::string& name, OptionsPtr opts) :
        Thresholds(Detector_t::Type_t::TAPSVeto,
                   BinSettings(100,0,20),
                   name, opts)
    {}
};
AUTO_REGISTER_PHYSICS(TAPSVeto_Thresholds)
