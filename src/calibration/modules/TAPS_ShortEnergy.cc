#include "TAPS_ShortEnergy.h"

#include "expconfig/detectors/TAPS.h"
#include "calibration/fitfunctions/FitLandau.h"
#include "calibration/fitfunctions/FitGaus.h"

#include "calibration/gui/CalCanvas.h"

#include "base/Logger.h"
#include "root-addons/cbtaps_display/TH2TAPS.h"

#include "TF1.h"


#include <list>
#include <cmath>

using namespace std;
using namespace ant;
using namespace ant::calibration;

TAPS_ShortEnergy::TAPS_ShortEnergy(
        const detector_ptr_t& taps,
        const std::shared_ptr<DataManager>& calmgr,
        Calibration::Converter::ptr_t converter,
        defaults_t defaultPedestals,
        defaults_t defaultGains,
        defaults_t defaultThresholds_Raw,
        defaults_t defaultThresholds_MeV,
        defaults_t defaultRelativeGains
        ) :
    Energy(taps,
           calmgr,
           converter,
           defaultPedestals,
           defaultGains,
           defaultThresholds_Raw,
           defaultThresholds_MeV,
           defaultRelativeGains,
           Channel_t::Type_t::IntegralShort),
    taps_detector(taps)
{

}

void TAPS_ShortEnergy::GetGUIs(std::list<std::unique_ptr<gui::CalibModule_traits> >& guis, OptionsPtr options)
{
    guis.emplace_back(std_ext::make_unique<GUI_Pedestals>(
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

TAPS_ShortEnergy::GUI_Gains::GUI_Gains(const string& basename,
                          OptionsPtr options,
                          CalibType& type,
                          const std::shared_ptr<DataManager>& calmgr,
                          const detector_ptr_t& taps) :
    GUI_CalibType(basename, options, type, calmgr, taps),
    func(make_shared<gui::FitGaus>()),
    taps_detector(taps)
{

}

TAPS_ShortEnergy::GUI_Gains::~GUI_Gains()
{

}

void TAPS_ShortEnergy::GUI_Gains::InitGUI(gui::ManagerWindow_traits& window)
{
    GUI_CalibType::InitGUI(window);

    canvas = window.AddCalCanvas();
    h_peaks = new TH1D("h_peaks","Peak positions",GetNumberOfChannels(),0,GetNumberOfChannels());
    h_peaks->SetXTitle("Channel Number");
    h_peaks->SetYTitle("Peak");
    h_relative = new TH1D("h_relative","Relative change from previous gains",GetNumberOfChannels(),0,GetNumberOfChannels());
    h_relative->SetXTitle("Channel Number");
    h_relative->SetYTitle("Relative change / %");

    h_relative_taps = new TH2TAPS("h_relative_taps",h_relative->GetTitle());
    h_peaks_taps = new TH2TAPS("h_peaks_taps",h_peaks->GetTitle());
}

std::shared_ptr<TH1> TAPS_ShortEnergy::GUI_Gains::GetHistogram(const WrapTFile& file) const
{
    return file.GetSharedHist<TH1>(CalibModule_traits::GetName()+"/rel_gamma");
}

gui::CalibModule_traits::DoFitReturn_t TAPS_ShortEnergy::GUI_Gains::DoFit(const TH1& hist, unsigned channel)
{
    if(detector->IsIgnored(channel) || taps_detector->IsPbWO4(channel))
        return DoFitReturn_t::Skip;

    auto& hist2 = dynamic_cast<const TH2&>(hist);

    h_projection = hist2.ProjectionX("h_projection",channel+1,channel+1);

    func->SetDefaults(h_projection);
    func->SetRange(interval<double>(-1,3));
    const auto it_fit_param = fitParameters.find(channel);
    if(it_fit_param != fitParameters.end()
       && !IgnorePreviousFitParameters) {
        VLOG(5) << "Loading previous fit parameters for channel " << channel;
        func->Load(it_fit_param->second);
    }

    for(size_t i=0;i<5;i++)
        func->Fit(h_projection);

    /// \todo implement automatic stop if fit failed?

    // goto next channel
    return DoFitReturn_t::Next;
}

void TAPS_ShortEnergy::GUI_Gains::DisplayFit()
{
    canvas->Divide(1,1);
    canvas->Show(h_projection, func.get());
}

void TAPS_ShortEnergy::GUI_Gains::StoreFit(unsigned channel)
{
    const double oldValue = previousValues[channel];

    const double convergenceFactor = 1.0;
    const double target_pos = 1.0;
    const double peak = func->GetPeakPosition();

    // apply convergenceFactor only to the desired procentual change of oldValue,
    // given by (pi0mass/pi0peak - 1)
    const double newValue = oldValue + oldValue * convergenceFactor * (target_pos/peak - 1);

    calibType.Values[channel] = newValue;

    const double relative_change = 100*(newValue/oldValue-1);

    LOG(INFO) << "Stored Ch=" << channel << ": PeakPosition " << peak
              << ", gain changed " << oldValue << " -> " << newValue
              << " (" << relative_change << " %)";


    // don't forget the fit parameters
    fitParameters[channel] = func->Save();

    h_peaks->SetBinContent(channel+1, peak);
    h_relative->SetBinContent(channel+1, relative_change);
}

bool TAPS_ShortEnergy::GUI_Gains::FinishSlice()
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


TAPS_ShortEnergy::GUI_Pedestals::GUI_Pedestals(const string& basename,
                                               OptionsPtr options,
                                               CalibType& type,
                                               const std::shared_ptr<DataManager>& calmgr,
                                               const detector_ptr_t& taps,
                                               std::shared_ptr<gui::PeakingFitFunction> fitfunction) :
    energy::GUI_Pedestals(basename, options, type, calmgr, taps, fitfunction),
    taps_detector(taps)
{
}

gui::CalibModule_traits::DoFitReturn_t TAPS_ShortEnergy::GUI_Pedestals::DoFit(const TH1& hist, unsigned channel)
{
    if(taps_detector->IsPbWO4(channel))
        return DoFitReturn_t::Skip;
    return energy::GUI_Pedestals::DoFit(hist, channel);
}
