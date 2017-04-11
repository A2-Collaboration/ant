#include "TAPS_Energy.h"

#include "Energy_GUI.h"

#include "expconfig/detectors/TAPS.h"

#include "calibration/fitfunctions/FitLandau.h"
#include "calibration/gui/CalCanvas.h"
#include "calibration/fitfunctions/FitGausPol3.h"

#include "root-addons/cbtaps_display/TH2TAPS.h"
#include "base/Logger.h"
#include "base/ParticleType.h"
#include "base/FloodFillAverages.h"

#include "TF1.h"

#include <list>
#include <cmath>

using namespace std;
using namespace ant;
using namespace ant::calibration;

vector<double> makeDefaults(const std::shared_ptr<const expconfig::detector::TAPS>& taps,
                         double default_BaF2,
                         double default_PbWO4)
{
    vector<double> defs(taps->GetNChannels(), default_BaF2);
    for(unsigned ch=0; ch<taps->GetNChannels(); ch++)
    {
        if(taps->IsPbWO4(ch))
        {
            defs[ch] = default_PbWO4;
        }
    }
    return defs;
}

TAPS_Energy::TAPS_Energy(const detector_ptr_t& taps,
        const std::shared_ptr<DataManager>& calmgr,
        const Calibration::Converter::ptr_t& converter,
        defaults_t defaultPedestals,
        defaults_t defaultGains,
        double defaultThreshold_Raw_BaF2,
        double defaultThreshold_Raw_PbWO4,
        defaults_t defaultThresholds_MeV,
        defaults_t defaultRelativeGains) :
    Energy(taps,
           calmgr,
           converter,
           defaultPedestals,
           defaultGains,
           makeDefaults(taps, defaultThreshold_Raw_BaF2, defaultThreshold_Raw_PbWO4),
           defaultThresholds_MeV,
           defaultRelativeGains),
    taps_detector(taps)
{
    // RelativeGains are flood filled
    RelativeGains.NotifyLoad = [taps] (CalibType& relativeGains) {
        auto& v = relativeGains.Values;
        auto getVal = [&v] (int ch) { return v[ch]; };
        auto setVal = [&v] (int ch, double val) {
            v[ch] = val;
            VLOG(5) << "Channel=" << ch << " flood filled";
        };
        auto getNeighbours = [taps] (int ch) { return taps->GetClusterElement(ch)->Neighbours; };
        auto getValid = [taps] (int ch) { return !taps->HasElementFlags(ch, Detector_t::ElementFlag_t::NoCalib); };
        floodFillAverages(v.size(), getVal, setVal, getNeighbours, getValid);
    };
}


void TAPS_Energy::GetGUIs(std::list<std::unique_ptr<gui::CalibModule_traits> >& guis, OptionsPtr options)
{
    guis.emplace_back(std_ext::make_unique<energy::GUI_Pedestals>(
                          GetName(),
                          options,
                          Pedestals,
                          calibrationManager,
                          taps_detector,
                          make_shared<gui::FitLandau>()
                          ));

    guis.emplace_back(std_ext::make_unique<GUI_Gains>(
                          GetName(),
                          options,
                          RelativeGains,
                          calibrationManager,
                          taps_detector
                          ));
}

struct FitTAPS_Energy : gui::FitGausPol3 {
    virtual void SetDefaults(TH1 *hist) override
    {
        // Amplitude
        func->SetParameter(0, 0.5*hist->GetMaximum());
        func->SetParLimits(0, 0, hist->GetMaximum());

        // x0
        auto range = GetRange();
        func->SetParameter(1, range.Clip(135.0));
        func->SetParLimits(1, range.Start(), range.Stop());

        // sigma
        func->SetParameter(2, 12);
        func->SetParLimits(2, 5, 50);

        func->SetParameter(3, 0.3*hist->GetMaximum());
        func->SetParameter(4, 0);
        func->SetParameter(5, 0);
        func->SetParameter(6, 0);

        Sync();
    }
};

TAPS_Energy::GUI_Gains::GUI_Gains(const string& basename, OptionsPtr options,
                          CalibType& type,
                          const std::shared_ptr<DataManager>& calmgr,
                          const std::shared_ptr<const expconfig::detector::TAPS>& taps_detector_) :
    GUI_CalibType(basename, options, type, calmgr, taps_detector_),
    func(make_shared<FitTAPS_Energy>()),
    taps_detector(taps_detector_)
{
}

void TAPS_Energy::GUI_Gains::InitGUI(gui::ManagerWindow_traits& window)
{
    GUI_CalibType::InitGUI(window);

    window.AddNumberEntry("Chi2/NDF limit for autostop", AutoStopOnChi2);
    window.AddNumberEntry("Minimum Fit Range", FitRange.Start());
    window.AddNumberEntry("Maximum Fit Range", FitRange.Stop());
    window.AddNumberEntry("Convergence Factor", ConvergenceFactor);

    canvas = window.AddCalCanvas();
    h_peaks = new TH1D("h_peaks","Peak positions",GetNumberOfChannels(),0,GetNumberOfChannels());
    h_peaks->SetXTitle("Channel Number");
    h_peaks->SetYTitle("Pi0 Peak / MeV");
    h_relative = new TH1D("h_relative","Relative change from previous gains",GetNumberOfChannels(),0,GetNumberOfChannels());
    h_relative->SetXTitle("Channel Number");
    h_relative->SetYTitle("Relative change / %");

    h_relative_taps = new TH2TAPS("h_relative_taps",h_relative->GetTitle());
    h_peaks_taps = new TH2TAPS("h_peaks_taps",h_peaks->GetTitle());
}

gui::CalibModule_traits::DoFitReturn_t TAPS_Energy::GUI_Gains::DoFit(const TH1& hist, unsigned channel)
{
    if(detector->IsIgnored(channel)) {
        VLOG(6) << "Skipping ignored channel " << channel;
        return DoFitReturn_t::Skip;
    }

    if(detector->HasElementFlags(channel, Detector_t::ElementFlag_t::NoCalib)) {
        VLOG(6) << "Skipping NoCalib-flagged channel " << channel;
        return DoFitReturn_t::Skip;
    }

    auto& hist2 = dynamic_cast<const TH2&>(hist);

    h_projection = hist2.ProjectionX("h_projection",channel+1,channel+1);

    // stop at empty histograms
    if(h_projection->GetEntries() < 1.0)
        return DoFitReturn_t::Display;

    func->SetDefaults(h_projection);
    func->SetRange(FitRange);
    const auto it_fit_param = fitParameters.find(channel);
    if(it_fit_param != fitParameters.end() && !IgnorePreviousFitParameters) {
        VLOG(5) << "Loading previous fit parameters for channel " << channel;
        func->Load(it_fit_param->second);
    }
    else {
        func->FitBackground(h_projection);
    }


    auto fit_loop = [this] (size_t retries) {

        const auto diff_at_side = .01;

        do {
            func->Fit(h_projection);
            VLOG(5) << "Chi2/dof = " << func->Chi2NDF();
            if(    (func->Chi2NDF() < AutoStopOnChi2)
                &&  func->EndsMatch(diff_at_side)
                ) {
                return true;
            }

             retries--;
        }
        while(retries>0);
        return false;
    };

    if(fit_loop(5))
        return DoFitReturn_t::Next;

    // try with defaults and background fit
    func->SetDefaults(h_projection);
    func->FitBackground(h_projection);

    if(fit_loop(5))
        return DoFitReturn_t::Next;

    // reached maximum retries without good chi2
    const auto range = func->GetRange();
    LOG(INFO) << "Chi2/dof = " << func->Chi2NDF() << " SBR_low = " << func->SignalToBackground(range.Start()) << " SBR_high = " << func->SignalToBackground(range.Stop());
    return DoFitReturn_t::Display;
}

void TAPS_Energy::GUI_Gains::DisplayFit()
{
    canvas->Divide(1,1);
    canvas->Show(h_projection, func.get());
}

void TAPS_Energy::GUI_Gains::StoreFit(unsigned channel)
{
    const double oldValue = previousValues[channel];
    const double pi0mass = ParticleTypeDatabase::Pi0.Mass();
    const double pi0peak = func->GetPeakPosition();

    // apply convergenceFactor only to the desired procentual change of oldValue,
    // given by (pi0mass/pi0peak - 1)
    const double newValue = oldValue + oldValue * ConvergenceFactor * (pi0mass/pi0peak - 1);

    calibType.Values[channel] = newValue;

    const double relative_change = 100*(newValue/oldValue-1);

    LOG(INFO) << "Stored Ch=" << channel << ": PeakPosition " << pi0peak
              << " MeV,  gain changed " << oldValue << " -> " << newValue
              << " (" << relative_change << " %)";


    // don't forget the fit parameters
    fitParameters[channel] = func->Save();

    h_peaks->SetBinContent(channel+1, pi0peak);
    h_relative->SetBinContent(channel+1, relative_change);
}

bool TAPS_Energy::GUI_Gains::FinishSlice()
{
    canvas->Clear();
    canvas->Divide(2,2);

    canvas->cd(1);
    h_peaks->SetStats(false);
    h_peaks->Draw("P");
    canvas->cd(2);
    h_peaks_taps->SetElements(*h_peaks);
    h_peaks_taps->Draw("colz");

    canvas->cd(3);
    h_relative->SetStats(false);
    h_relative->Draw("P");
    canvas->cd(4);
    h_relative_taps->SetElements(*h_relative);
    h_relative_taps->Draw("colz");

    return true;
}
