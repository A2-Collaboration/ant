#include "Energy_GUI.h"

#include "calibration/DataManager.h"
#include "calibration/gui/CalCanvas.h"
#include "calibration/fitfunctions/FitLandauPol0.h"
#include "calibration/fitfunctions/FitLandauExpo.h"
#include "calibration/fitfunctions/FitWeibullLandauPol1.h"
#include "calibration/fitfunctions/FitVetoBand.h"

#include "tree/TCalibrationData.h"

#include "base/math_functions/Linear.h"
#include "base/std_ext/math.h"
#include "base/std_ext/string.h"
#include "base/TH_ext.h"
#include "base/Logger.h"

#include "TH2.h"
#include "TH3.h"
#include "TF1.h"

#include "TGNumberEntry.h"

using namespace std;
using namespace ant;
using namespace ant::calibration;
using namespace ant::calibration::energy;

GUI_Pedestals::GUI_Pedestals(const string& basename,
        OptionsPtr options,
        CalibType& type,
        const std::shared_ptr<DataManager>& calmgr,
        const detector_ptr_t& detector,
        shared_ptr<gui::PeakingFitFunction> fitfunction) :
    GUI_CalibType(basename, options, type, calmgr, detector, Calibration::AddMode_t::RightOpen),
    func(fitfunction)
{

}

void GUI_Pedestals::InitGUI(gui::ManagerWindow_traits& window)
{
    GUI_CalibType::InitGUI(window);
    canvas = window.AddCalCanvas();
}

gui::CalibModule_traits::DoFitReturn_t GUI_Pedestals::DoFit(const TH1& hist, unsigned channel)
{
    if(detector->IsIgnored(channel))
        return DoFitReturn_t::Skip;

    auto& hist2 = dynamic_cast<const TH2&>(hist);

    h_projection = hist2.ProjectionX("h_projection",channel+1,channel+1);

    func->SetDefaults(h_projection);
    const auto it_fit_param = fitParameters.find(channel);
    if(it_fit_param != fitParameters.end() && !IgnorePreviousFitParameters) {
        VLOG(5) << "Loading previous fit parameters for channel " << channel;
        func->Load(it_fit_param->second);
    }

    for(size_t i=0;i<5;i++)
        func->Fit(h_projection);

    /// \todo implement automatic stop if fit failed?

    // goto next channel
    return DoFitReturn_t::Next;
}

void GUI_Pedestals::DisplayFit()
{
    canvas->Show(h_projection, func.get());
}

void GUI_Pedestals::StoreFit(unsigned channel)
{

    const double oldValue = previousValues[channel];
    const double newValue = func->GetPeakPosition();

    calibType.Values[channel] = newValue;

    const double relative_change = 100*(newValue/oldValue-1);

    LOG(INFO) << "Stored Ch=" << channel << ":  "
              <<" Pedestal changed " << oldValue << " -> " << newValue
              << " (" << relative_change << " %)";


    // don't forget the fit parameters
    fitParameters[channel] = func->Save();
}

bool GUI_Pedestals::FinishSlice()
{
    return false;
}

GUI_Banana::GUI_Banana(const string& basename,
                               OptionsPtr options,
                               CalibType& type,
                               const std::shared_ptr<DataManager>& calmgr,
                               const detector_ptr_t& detector,
                               const interval<double>& projectionrange,
                               const double proton_peak_mc_pos
                               ) :
    GUI_CalibType(basename, options, type, calmgr, detector),
    func(make_shared<gui::FitLandauPol0>()),
    projection_range(projectionrange),
    proton_peak_mc(proton_peak_mc_pos),
    full_hist_name(
            options->Get<string>("HistogramPath", CalibModule_traits::GetName())
            + "/"
            + options->Get<string>("HistogramName", "Bananas"))
{

}

std::shared_ptr<TH1> GUI_Banana::GetHistogram(const WrapTFile& file) const
{
    return file.GetSharedHist<TH1>(full_hist_name);
}

void GUI_Banana::InitGUI(gui::ManagerWindow_traits& window)
{
    GUI_CalibType::InitGUI(window);
    window.AddNumberEntry("Chi2/NDF limit for autostop", AutoStopOnChi2);

    c_fit = window.AddCalCanvas();
    c_extra = window.AddCalCanvas();


    h_relative = new TH1D("h_relative","Relative change from previous gains",GetNumberOfChannels(),0,GetNumberOfChannels());
    h_relative->SetXTitle("Channel Number");
    h_relative->SetYTitle("Relative change / %");
}

gui::CalibModule_traits::DoFitReturn_t GUI_Banana::DoFit(const TH1& hist, unsigned ch)
{
    if(detector->IsIgnored(ch))
        return DoFitReturn_t::Skip;

    auto& h_bananas = dynamic_cast<const TH3&>(hist);
    h_bananas.GetZaxis()->SetRange(ch+1,ch+1);
    banana = dynamic_cast<TH2D*>(h_bananas.Project3D("yx"));
    auto xaxis = banana->GetXaxis();
    h_projection = dynamic_cast<TH1D*>(banana->ProjectionY(
                                           "_py",
                                           xaxis->FindFixBin(projection_range.Start()),
                                           xaxis->FindFixBin(projection_range.Stop())
                                           )
                                       );

    // stop at empty histograms
    if(h_projection->GetEntries()==0)
        return DoFitReturn_t::Display;

    func->SetRange(interval<double>(0.3,6));
    func->SetDefaults(h_projection);
    const auto it_fit_param = fitParameters.find(ch);
    if(it_fit_param != fitParameters.end() && !IgnorePreviousFitParameters) {
        VLOG(5) << "Loading previous fit parameters for channel " << ch;
        func->Load(it_fit_param->second);
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

    // reached maximum retries without good chi2
    LOG(INFO) << "Chi2/dof = " << func->Chi2NDF();
    return DoFitReturn_t::Display;
}

void GUI_Banana::DisplayFit()
{
    c_fit->Show(h_projection, func.get());

    c_extra->cd();
    banana->Draw("colz");
}

void GUI_Banana::StoreFit(unsigned channel)
{
    const double oldValue = previousValues[channel];

    const double protonpeak = func->GetPeakPosition();

    const double newValue = oldValue * proton_peak_mc / protonpeak;

    calibType.Values[channel] = newValue;

    const double relative_change = 100*(newValue/oldValue-1);

    LOG(INFO) << "Stored Ch=" << channel << ": Peak " << protonpeak
              << " MeV,  gain changed " << oldValue << " -> " << newValue
              << " (" << relative_change << " %)";


    // don't forget the fit parameters
    fitParameters[channel] = func->Save();

    h_relative->SetBinContent(channel+1, relative_change);

}

bool GUI_Banana::FinishSlice()
{
    c_extra->Clear();
    c_fit->Clear();

    c_fit->cd();
    h_relative->SetStats(false);
    h_relative->Draw("P");

    return true;
}


GUI_HEP::GUI_HEP(const string& basename,
                         OptionsPtr options,
                         CalibType& type,
                         const std::shared_ptr<DataManager>& calmgr,
                         const detector_ptr_t& detector,
                         const double proton_peak_mc_pos
                         ) :
    GUI_CalibType(basename, options, type, calmgr, detector),
    func(make_shared<gui::FitWeibullLandauPol1>()),
    proton_peak_mc(proton_peak_mc_pos),
    full_hist_name(
            options->Get<string>("HistogramPath", CalibModule_traits::GetName())
            + "/"
            + options->Get<string>("HistogramName", "projections_hep"))
{

}

std::shared_ptr<TH1> GUI_HEP::GetHistogram(const WrapTFile& file) const
{
    return file.GetSharedHist<TH1>(full_hist_name);
}

void GUI_HEP::InitGUI(gui::ManagerWindow_traits& window)
{
    GUI_CalibType::InitGUI(window);
    window.AddNumberEntry("Chi2/NDF limit for autostop", AutoStopOnChi2);

    canvas = window.AddCalCanvas();

    h_peaks = new TH1D("h_peaks","Peak positions",GetNumberOfChannels(),0,GetNumberOfChannels());
    h_peaks->SetXTitle("Channel Number");
    h_peaks->SetYTitle("High Energy Proton Peak / MeV");
    h_relative = new TH1D("h_relative","Relative change from previous gains",GetNumberOfChannels(),0,GetNumberOfChannels());
    h_relative->SetXTitle("Channel Number");
    h_relative->SetYTitle("Relative change / %");

    LOG(INFO) << "Use high energy protons for PID gain calibration";
}

gui::CalibModule_traits::DoFitReturn_t GUI_HEP::DoFit(const TH1& hist, unsigned ch)
{
    if(detector->IsIgnored(ch))
        return DoFitReturn_t::Skip;

    auto& hist2 = dynamic_cast<const TH2&>(hist);
    h_projection = hist2.ProjectionX("h_projection",ch+1,ch+1);

    // stop at empty histograms
    if(h_projection->GetEntries()==0)
        return DoFitReturn_t::Display;

    //auto range = interval<double>(1,9);
    auto range = interval<double>(.8,9.5);

    func->SetDefaults(h_projection);
    func->SetRange(range);
    const auto it_fit_param = fitParameters.find(ch);
    if(it_fit_param != fitParameters.end() && !IgnorePreviousFitParameters && false) {
        VLOG(5) << "Loading previous fit parameters for channel " << ch;
        func->Load(it_fit_param->second);
    }
    else {
        func->FitSignal(h_projection);
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

    // try with defaults ...
    func->SetDefaults(h_projection);
    func->Fit(h_projection);

    if(fit_loop(5))
        return DoFitReturn_t::Next;

    // ... and with defaults and first a signal only fit ...
    func->SetDefaults(h_projection);
    func->FitSignal(h_projection);

    if(fit_loop(5))
        return DoFitReturn_t::Next;

    // ... and as a last resort background, signal and a few last fit tries
    func->SetDefaults(h_projection);
    func->FitBackground(h_projection);
    func->Fit(h_projection);
    func->FitSignal(h_projection);

    if(fit_loop(5))
        return DoFitReturn_t::Next;

    // reached maximum retries without good chi2
    LOG(INFO) << "Chi2/dof = " << func->Chi2NDF();
    return DoFitReturn_t::Display;
}

void GUI_HEP::DisplayFit()
{
    canvas->Show(h_projection, func.get());
}

void GUI_HEP::StoreFit(unsigned channel)
{
    const double oldValue = previousValues[channel];
    const double peak = func->GetPeakPosition();
    const double newValue = oldValue * proton_peak_mc / peak;

    calibType.Values[channel] = newValue;

    const double relative_change = 100*(newValue/oldValue-1);

    LOG(INFO) << "Stored Ch=" << channel << ": PeakPosition " << peak
              << " MeV,  gain changed " << oldValue << " -> " << newValue
              << " (" << relative_change << " %)";


    // don't forget the fit parameters
    fitParameters[channel] = func->Save();

    h_peaks->SetBinContent(channel+1, peak);
    h_relative->SetBinContent(channel+1, relative_change);

}

bool GUI_HEP::FinishSlice()
{
    canvas->Clear();
    canvas->Divide(1,2);

    canvas->cd(1);
    h_peaks->SetStats(false);
    h_peaks->Draw("P");
    canvas->cd(2);
    h_relative->SetStats(false);
    h_relative->Draw("P");

    return true;
}


GUI_BananaSlices::GUI_BananaSlices(const string& basename,
                         OptionsPtr options,
                         CalibType& type,
                         const std::shared_ptr<DataManager>& calmgr,
                         const detector_ptr_t& detector,
                         const interval<double>& fitrange
                         ) :
    GUI_CalibType(basename, options, type, calmgr, detector),
    func(make_shared<gui::FitVetoBand>()),
    fit_range(fitrange),
    full_hist_name(
        options->Get<string>("HistogramPath", CalibModule_traits::GetName())
        + "/"
        + options->Get<string>("HistogramName", "Bananas"))
{
    slicesY_gaus = new TF1("slicesY_gaus","gaus");
}

std::shared_ptr<TH1> GUI_BananaSlices::GetHistogram(const WrapTFile& file) const
{
    return file.GetSharedHist<TH1>(full_hist_name);
}

void GUI_BananaSlices::InitGUI(gui::ManagerWindow_traits& window)
{
    GUI_CalibType::InitGUI(window);

    c_fit = window.AddCalCanvas();
    c_extra = window.AddCalCanvas();

    window.AddNumberEntry("Lower Energy Limit for fit function",fit_range.Start(),[this] (const TGNumberEntry& e) {
        fit_range.Start() = e.GetNumber();
        func->SetRange(fit_range);
        c_fit->UpdateMe();
    });

    window.AddNumberEntry("Upper Energy Limit for fit function",fit_range.Stop(),[this] (const TGNumberEntry& e) {
        fit_range.Stop() = e.GetNumber();
        func->SetRange(fit_range);
        c_fit->UpdateMe();
    });

    window.AddNumberEntry("Chi2/NDF limit for autostop", AutoStopOnChi2);
    window.AddNumberEntry("SlicesYEntryCut", slicesY_entryCut);
    window.AddNumberEntry("SlicesYIQRFactor low  (outlier detection)", slicesY_IQRFactor_lo);
    window.AddNumberEntry("SlicesYIQRFactor high (outlier detection)", slicesY_IQRFactor_hi);

    h_vals = new TH1D("h_vals","Energy values from Veto band",GetNumberOfChannels(),0,GetNumberOfChannels());
    h_vals->SetXTitle("Channel Number");
    h_vals->SetYTitle("Calculated Veto Energy / MeV");
    h_relative = new TH1D("h_relative","Relative change from previous gains",GetNumberOfChannels(),0,GetNumberOfChannels());
    h_relative->SetXTitle("Channel Number");
    h_relative->SetYTitle("Relative change / %");

    LOG(INFO) << "Use proton bananas for PID gain calibration";
    LOG(WARNING) << "Please make sure to set a fixed energy fitting range"
                 << " via the GUI number fields and keep it for all channels!";
}

gui::CalibModule_traits::DoFitReturn_t GUI_BananaSlices::DoFit(const TH1& hist, unsigned ch)
{
    if(detector->IsIgnored(ch))
        return DoFitReturn_t::Skip;

    auto& h_vetoband = dynamic_cast<const TH3&>(hist);

    h_vetoband.GetZaxis()->SetRange(ch+1,ch+1);
    h_proj = dynamic_cast<TH2D*>(h_vetoband.Project3D("yx"));

    h_means = TH_ext::FitSlicesY(h_proj, slicesY_gaus, slicesY_entryCut,
                                 slicesY_IQRFactor_lo, slicesY_IQRFactor_hi);
    //h_means->SetMinimum(proj.GetYaxis()->GetXmin());
    //h_means->SetMaximum(proj.GetYaxis()->GetXmax());

    // stop at empty histograms
    if(h_means->GetEntries()==0)
        return DoFitReturn_t::Display;

    func->SetDefaults(h_means);
    func->SetRange(fit_range);
    const auto it_fit_param = fitParameters.find(ch);
    if(it_fit_param != fitParameters.end()) {
        VLOG(5) << "Loading previous fit parameters for channel " << ch;
        func->Load(it_fit_param->second);
        func->SetRange(fit_range);
    }


    auto fit_loop = [this] (size_t retries) {
        do {
            func->Fit(h_means);
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

    // if the fit failed, use defaults and adjust offset by fitting background, let the user handle the rest
    func->SetDefaults(h_means);
    func->FitBackground(h_means);

    // reached maximum retries without good chi2
    LOG(INFO) << "Chi2/dof = " << func->Chi2NDF();
    return DoFitReturn_t::Display;
}

void GUI_BananaSlices::DisplayFit()
{
    c_fit->Show(h_means, func.get(), true);

    c_extra->cd();
    h_proj->DrawCopy("colz");
    func->Draw();
}

void GUI_BananaSlices::StoreFit(unsigned channel)
{
    const double energy = fit_range.Stop();
    const double oldValue = previousValues[channel];
    const double val = func->Eval(energy);
    const double ref = func->EvalReference(energy);
    const double newValue = oldValue * ref/val;

    calibType.Values[channel] = newValue;

    const double relative_change = 100*(newValue/oldValue-1);

    LOG(INFO) << "Stored Ch=" << channel << ": Energy value at " << energy
              << " MeV: " << val << " MeV, reference: " << ref << " MeV"
              << " ;  gain changed " << oldValue << " -> " << newValue
              << " (" << relative_change << " %)";


    // don't forget the fit parameters
    fitParameters[channel] = func->Save();

    h_vals->SetBinContent(channel+1, val);
    h_relative->SetBinContent(channel+1, relative_change);

    //LOG(INFO) << "Stored Ch=" << channel << " Parameters: " << fitParameters[channel];
}

bool GUI_BananaSlices::FinishSlice()
{
    // don't request stop...
    return false;

//    canvas->Clear();
//    canvas->Divide(1,2);

//    canvas->cd(1);
//    h_vals->SetStats(false);
//    h_vals->Draw("P");
//    canvas->cd(2);
//    h_relative->SetStats(false);
//    h_relative->Draw("P");

//    return true;
}
