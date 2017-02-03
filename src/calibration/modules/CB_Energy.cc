#include "CB_Energy.h"

#include "calibration/gui/CalCanvas.h"
#include "calibration/fitfunctions/FitGausPol3.h"

#include "expconfig/detectors/CB.h"

#include "root-addons/cbtaps_display/TH2CB.h"
#include "base/Logger.h"
#include "base/ParticleType.h"

#include <list>


using namespace std;
using namespace ant;
using namespace ant::calibration;

CB_Energy::CB_Energy(std::shared_ptr<expconfig::detector::CB> cb,
                     std::shared_ptr<DataManager> calmgr,
                     Calibration::Converter::ptr_t converter,
                     const std::vector<double>& defaultPedestals,
                     const std::vector<double>& defaultGains,
                     const std::vector<double>& defaultThresholds,
                     const std::vector<double>& defaultRelativeGains):
    Energy(cb->Type,
           calmgr,
           converter,
           defaultPedestals,
           defaultGains,
           defaultThresholds,
           defaultRelativeGains),
    cb_detector(cb)
{

}

void CB_Energy::GetGUIs(list<unique_ptr<gui::CalibModule_traits> >& guis, OptionsPtr options)
{
    guis.emplace_back(std_ext::make_unique<GUI_Gains>(
                          GetName(),
                          options,
                          RelativeGains,
                          calibrationManager,
                          cb_detector
                          ));
}

CB_Energy::GUI_Gains::GUI_Gains(const string& basename,
                          OptionsPtr options,
                          CalibType& type,
                          const std::shared_ptr<DataManager>& calmgr,
                          const std::shared_ptr<Detector_t>& detector) :
    GUI_CalibType(basename, options, type, calmgr, detector),
    func(make_shared<gui::FitGausPol3>())
{
}

void CB_Energy::GUI_Gains::InitGUI(gui::ManagerWindow_traits* window)
{
    GUI_CalibType::InitGUI(window);

    window->AddNumberEntry("Chi2/NDF limit for autostop", AutoStopOnChi2);
    window->AddNumberEntry("Minimum Fit Range", FitRange.Start());
    window->AddNumberEntry("Maximum Fit Range", FitRange.Stop());
    window->AddNumberEntry("Convergence Factor", ConvergenceFactor);

    canvas = window->AddCalCanvas();
    h_peaks = new TH1D("h_peaks","Peak positions",GetNumberOfChannels(),0,GetNumberOfChannels());
    h_peaks->SetXTitle("Channel Number");
    h_peaks->SetYTitle("Pi0 Peak / MeV");
    h_relative = new TH1D("h_relative","Relative change from previous gains",GetNumberOfChannels(),0,GetNumberOfChannels());
    h_relative->SetXTitle("Channel Number");
    h_relative->SetYTitle("Relative change / %");

    h_relative_cb = new TH2CB("h_relative_cb",h_relative->GetTitle());
    h_peaks_cb = new TH2CB("h_peaks_cb",h_peaks->GetTitle());
}

gui::CalibModule_traits::DoFitReturn_t CB_Energy::GUI_Gains::DoFit(TH1* hist, unsigned channel)
{
    if(detector->IsIgnored(channel))
        return DoFitReturn_t::Skip;

    TH2* hist2 = dynamic_cast<TH2*>(hist);

    h_projection = hist2->ProjectionX("h_projection",channel+1,channel+1);

    // stop at empty histograms
    if(h_projection->GetEntries()==0)
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
        do {
            func->Fit(h_projection);
            VLOG(5) << "Chi2/dof = " << func->Chi2NDF();
            if(func->Chi2NDF() < AutoStopOnChi2) {
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
    LOG(INFO) << "Chi2/dof = " << func->Chi2NDF();
    return DoFitReturn_t::Display;
}

void CB_Energy::GUI_Gains::DisplayFit()
{
    canvas->Divide(1,1);
    canvas->Show(h_projection, func.get());
}

void CB_Energy::GUI_Gains::StoreFit(unsigned channel)
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

bool CB_Energy::GUI_Gains::FinishSlice()
{
    canvas->Clear();
    canvas->Divide(2,2);

    canvas->cd(1);
    h_peaks->SetStats(false);
    h_peaks->Draw("P");
    canvas->cd(2);
    h_peaks_cb->SetElements(*h_peaks);
    h_peaks_cb->Draw("colz");

    canvas->cd(3);
    h_relative->SetStats(false);
    h_relative->Draw("P");
    canvas->cd(4);
    h_relative_cb->SetElements(*h_relative);
    h_relative_cb->Draw("colz");

    return true;
}

