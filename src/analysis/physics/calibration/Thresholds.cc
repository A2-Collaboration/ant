#include "Thresholds.h"

#include "expconfig/ExpConfig.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

Thresholds::Thresholds(const Detector_t::Type_t& detectorType,
                       const BinSettings& bins_x,
                       const string& name, OptionsPtr opts) :
    Physics(name, opts)
{
    auto setup = ExpConfig::Setup::GetLastFound();
    Detector = setup->GetDetector(detectorType);

    hThresholds_ADC = HistFac.makeTH2D("Thresholds ADC","Energy","Element",
                                       bins_x,
                                       BinSettings(Detector->GetNChannels()),
                                       "hThresholds_ADC");

    hThresholds_TDC = HistFac.makeTH2D("Thresholds TDC","Energy","Element",
                                       bins_x,
                                       BinSettings(Detector->GetNChannels()),
                                       "hThresholds_TDC");
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
    }
}

void Thresholds::ShowResult()
{
    canvas(GetName()) << drawoption("colz")
                      << hThresholds_ADC
                      << hThresholds_TDC
                      << endc;
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
