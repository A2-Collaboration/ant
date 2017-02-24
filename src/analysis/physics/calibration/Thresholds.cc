#include "Thresholds.h"

#include "expconfig/ExpConfig.h"
#include "root-addons/cbtaps_display/TH2CB.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

Thresholds::Thresholds(const Detector_t::Type_t& detectorType,
                       const BinSettings& bins_x,
                       const string& name, OptionsPtr opts) :
    Physics(name, opts),
    Detector(ExpConfig::Setup::GetDetector(detectorType))
{
    const BinSettings bins_ch(Detector->GetNChannels());

    hThresholds_ADC = HistFac.makeTH2D("Thresholds ADC","Energy","Element",
                                       bins_x, bins_ch,
                                       "hThresholds_ADC");

    hThresholds_TDC = HistFac.makeTH2D("Thresholds TDC","Energy","Element",
                                       bins_x, bins_ch,
                                       "hThresholds_TDC");
    // makes only sense if Ant is running with '-S IncludeIgnoredElements=1'
    hMaybeDeadTDCs = HistFac.makeTH1D("ADC>4.5MeV but no TDC","Channel","",
                                      bins_ch, "hMaybeDeadTDCs");
}

void Thresholds::ProcessEvent(const TEvent& event, manager_t&)
{
    struct hit_t {
        double Energy = std_ext::NaN;
        double Time = std_ext::NaN;
    };

    map<unsigned, hit_t> hits;
    for(const TDetectorReadHit& readhit : event.Reconstructed().DetectorReadHits) {
        if(readhit.DetectorType != Detector->Type)
            continue;
        if(readhit.Values.empty())
            continue;
        auto& hit = hits[readhit.Channel];
        if(readhit.ChannelType == Channel_t::Type_t::Integral)
            hit.Energy = readhit.Values.front();
        if(readhit.ChannelType == Channel_t::Type_t::Timing)
            hit.Time = readhit.Values.front();
    }
    for(const auto& it_hit : hits) {
        auto ch = it_hit.first;
        const hit_t& hit = it_hit.second;
        hThresholds_ADC->Fill(hit.Energy, ch);
        if(isfinite(hit.Time))
            hThresholds_TDC->Fill(hit.Energy, ch);
        else if(hit.Energy > 4.5)
            hMaybeDeadTDCs->Fill(ch);
    }
}

void Thresholds::ShowResult()
{
    canvas c(GetName());
    c << drawoption("colz")
      << hThresholds_ADC
      << hThresholds_TDC
      << hMaybeDeadTDCs;
    if(Detector->Type == Detector_t::Type_t::CB) {
        auto cbhist = new TH2CB("cbhist", hMaybeDeadTDCs->GetTitle());
        cbhist->FillElements(*hMaybeDeadTDCs);
        for(unsigned ch=0;ch<Detector->GetNChannels();ch++) {
            if(Detector->IsIgnored(ch)) {
                cbhist->CreateMarker(ch);
            }
        }
        c << cbhist;
    }
    c << endc;
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

struct TAPS_Thresholds : Thresholds {
    TAPS_Thresholds(const std::string& name, OptionsPtr opts) :
        Thresholds(Detector_t::Type_t::TAPS,
                   BinSettings(100,0,20),
                   name, opts)
    {}
};
AUTO_REGISTER_PHYSICS(TAPS_Thresholds)
