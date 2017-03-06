#include "Energy.h"

#include "calibration/DataManager.h"
#include "calibration/gui/CalCanvas.h"
#include "calibration/fitfunctions/FitGausPol1.h"
#include "calibration/fitfunctions/FitLandauExpo.h"
#include "calibration/fitfunctions/FitWeibullLandauPol1.h"

#include "tree/TCalibrationData.h"
#include "tree/TDetectorReadHit.h"

#include "base/Logger.h"
#include "base/std_ext/math.h"
#include "base/math_functions/Linear.h"

#include "TH2.h"
#include "TH3.h"
#include "TF1.h"

#include <cstdint>
#include <ctime>
#include <algorithm>
#include <sstream>
#include <vector>
#include <list>
#include <cmath>

using namespace std;
using namespace ant;
using namespace ant::calibration;
using namespace ant::std_ext;

Energy::Energy(const detector_ptr_t& det,
               const std::shared_ptr<DataManager>& calmgr,
               const Calibration::Converter::ptr_t& converter,
               defaults_t defaultPedestals,
               defaults_t defaultGains,
               defaults_t defaultThresholds_Raw,
               defaults_t defaultThresholds_MeV,
               defaults_t defaultRelativeGains,
               Channel_t::Type_t channelType) :
    Calibration::Module(
        std_ext::formatter()
        << Detector_t::ToString(det->Type)
        << "_"
        << ( channelType == Channel_t::Type_t::IntegralShort ? "Short" : "" )
        << "Energy"
           ),
    DetectorType(det->Type),
    ChannelType(channelType),
    calibrationManager(calmgr),
    Converter(move(converter)),
    Pedestals(det, "Pedestals", defaultPedestals),
    Gains(det, "Gains", defaultGains, "ggIM"),
    Thresholds_Raw(det, "Thresholds_Raw", defaultThresholds_Raw),
    Thresholds_MeV(det, "Thresholds_MeV", defaultThresholds_MeV),
    RelativeGains(det, "RelativeGains", defaultRelativeGains, "ggIM")
{
    if(Converter==nullptr)
        throw std::runtime_error("Given converter should not be nullptr");
}

Energy::~Energy()
{
}

void Energy::ApplyTo(const readhits_t& hits)
{
    const auto& dethits = hits.get_item(DetectorType);

    // now calibrate the Energies (ignore any other kind of hits)
    for(TDetectorReadHit& dethit : dethits) {
        if(dethit.ChannelType != ChannelType)
            continue;

        // prefer building from RawData if available
        if(!dethit.RawData.empty()) {
            // clear previously read values (if any)
            dethit.Values.resize(0);

            // apply pedestal/gain to each of the values (might be multihit)
            for(const double& conv : Converter->Convert(dethit.RawData)) {
                TDetectorReadHit::Value_t value(conv);
                value.Calibrated -= Pedestals.Get(dethit.Channel);

                const double threshold = Thresholds_Raw.Get(dethit.Channel);
                if(value.Calibrated<threshold)
                    continue;

                // calibrate with absolute gain
                value.Calibrated *= Gains.Get(dethit.Channel);

                dethit.Values.emplace_back(move(value));
            }
        }

        // apply relative gain and threshold on MC
        {
            auto it_value = dethit.Values.begin();
            while(it_value != dethit.Values.end()) {
                it_value->Calibrated *= RelativeGains.Get(dethit.Channel);

                if(IsMC) {
                    const double threshold = Thresholds_MeV.Get(dethit.Channel);
                    // erase from Values if below threshold
                    if(it_value->Calibrated<threshold) {
                        it_value = dethit.Values.erase(it_value);
                        continue;
                    }
                }

                ++it_value;
            }
        }
    }
}

double Energy::CalibType::Get(unsigned channel) const {

    if(Values.empty()) {
        if(DefaultValues.size() == 1) {
            return DefaultValues.front();
        }
        else {
            return DefaultValues.at(channel);
        }
    }
    else {
        return Values.at(channel);
    }
}

Energy::CalibType::CalibType(
        const std::shared_ptr<const Detector_t>& det,
        const string& name,
        const std::vector<double>& defaultValues,
        const string& histname) :
    Name(name),
    // use name for histogram if not provided different
    HistogramName(histname.empty() ? name : histname),
    Values(),
    DefaultValues(defaultValues)
{
    if(DefaultValues.size() != 1 && DefaultValues.size() != det->GetNChannels()) {
        throw runtime_error("Wrong size of default values for calibType="+name+" det="+Detector_t::ToString(det->Type));
    }
}

std::list<Updateable_traits::Loader_t> Energy::GetLoaders()
{

    std::list<Updateable_traits::Loader_t> loaders;

    for(auto calibration : AllCalibrations) {

        auto loader = [this, calibration]
                (const TID& currPoint, TID& nextChangePoint)
        {
            TCalibrationData cdata;
            if(calibrationManager->GetData(
                   GetName()+"_"+ calibration->Name,
                   currPoint, cdata, nextChangePoint))
            {
                auto& values = calibration->Values;
                for (const auto& val: cdata.Data) {
                    if(values.size()<val.Key+1)
                        values.resize(val.Key+1);
                    values[val.Key] = val.Value;
                }

                // call notify load if present
                if(calibration->NotifyLoad)
                    calibration->NotifyLoad(*calibration);
            }
            else {
                LOG_IF(!calibration->Values.empty(), WARNING)
                        << "No calibration data found for " << calibration->Name
                        << " at changepoint TID=" << currPoint << ", using default values";
                calibration->Values.resize(0);
            }
        };

        loaders.emplace_back(loader);
    }

    return loaders;
}

void Energy::UpdatedTIDFlags(const TID& id)
{
    IsMC = id.isSet(TID::Flags_t::MC);
}

Energy::GUI_CalibType::GUI_CalibType(const string& basename, OptionsPtr opts,
                                     CalibType& type,
                                     const shared_ptr<DataManager>& calmgr,
                                     const detector_ptr_t& detector_,
                                     Calibration::AddMode_t mode) :
    gui::CalibModule_traits(basename),
    options(opts),
    calibType(type),
    calibrationManager(calmgr),
    detector(detector_),
    addMode(mode)
{}

string Energy::GUI_CalibType::GetName() const
{
    // serves as the CalibrationID for the manager
    return  CalibModule_traits::GetName()+"_"+calibType.Name;
}

shared_ptr<TH1> Energy::GUI_CalibType::GetHistogram(const WrapTFile& file) const
{
    // histogram name created by the specified Physics class
    return file.GetSharedHist<TH1>(options->Get<string>("HistogramPath", CalibModule_traits::GetName()) + "/"+calibType.HistogramName);
}

unsigned Energy::GUI_CalibType::GetNumberOfChannels() const
{
    return detector->GetNChannels();
}

void Energy::GUI_CalibType::InitGUI(gui::ManagerWindow_traits* window) {
    window->AddCheckBox("Ignore prev fit params", IgnorePreviousFitParameters);
    window->AddCheckBox("Use params from prev slice", UsePreviousSliceParams);
}

void Energy::GUI_CalibType::StartSlice(const interval<TID>& range)
{
    // always make sure the values are large enough
    std::vector<double>  values;
    values.resize(GetNumberOfChannels());
    for(size_t i=0; i<values.size(); ++i) {
        values[i] = calibType.Get(i);
    }

    TCalibrationData cdata;
    if(calibrationManager->GetData(GetName(), range.Start(), cdata)) {
        for(const TKeyValue<double>& kv : cdata.Data) {
            values[kv.Key] = kv.Value;
        }
        LOG(INFO) << GetName() << ": Loaded previous values from database";

        if(fitParameters.empty() || !UsePreviousSliceParams) {
            for(const TKeyValue<vector<double>>& kv : cdata.FitParameters) {
                fitParameters[kv.Key] = kv.Value;
            }
            LOG(INFO) << GetName() << ": Loaded previous fit parameter from database";
        }
        else if(!fitParameters.empty()) {
            LOG(INFO) << GetName() << ": Using fit parameters from previous slice";
        }
    }
    else {
        LOG(INFO) << GetName() << ": No previous values found, built from default value";
    }

    calibType.Values = values;

    // save a copy for comparison at finish stage
    previousValues = calibType.Values;

}

void Energy::GUI_CalibType::StoreFinishSlice(const interval<TID>& range)
{
    TCalibrationData cdata(
                GetName(),
                range.Start(),
                range.Stop()
                );

    std::vector<double>& values = calibType.Values;

    // fill data
    for(unsigned ch=0;ch<values.size();ch++) {
        cdata.Data.emplace_back(ch, values[ch]);
    }

    // fill fit parameters (if any)
    for(const auto& it_map : fitParameters) {
        const unsigned ch = it_map.first;
        const vector<double>& params = it_map.second;
        cdata.FitParameters.emplace_back(ch, params);
    }

    calibrationManager->Add(cdata, addMode);
}

Energy::GUI_Pedestals::GUI_Pedestals(const string& basename,
        OptionsPtr options,
        CalibType& type,
        const std::shared_ptr<DataManager>& calmgr,
        const detector_ptr_t& detector,
        shared_ptr<gui::PeakingFitFunction> fitfunction) :
    GUI_CalibType(basename, options, type, calmgr, detector, Calibration::AddMode_t::RightOpen),
    func(fitfunction)
{

}

void Energy::GUI_Pedestals::InitGUI(gui::ManagerWindow_traits* window)
{
    GUI_CalibType::InitGUI(window);
    canvas = window->AddCalCanvas();
}

gui::CalibModule_traits::DoFitReturn_t Energy::GUI_Pedestals::DoFit(const TH1& hist, unsigned channel)
{
    if(detector->IsIgnored(channel))
        return DoFitReturn_t::Skip;

    auto hist2 = dynamic_cast<const TH2&>(hist);

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

void Energy::GUI_Pedestals::DisplayFit()
{
    canvas->Show(h_projection, func.get());
}

void Energy::GUI_Pedestals::StoreFit(unsigned channel)
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

bool Energy::GUI_Pedestals::FinishSlice()
{
    return false;
}

struct FitProtonPeak : gui::FitGausPol1 {
    virtual void SetDefaults(TH1 *hist) override {

        const auto range = GetRange();

        const auto startbin = hist->FindBin(range.Start());
        const auto stopbin  = hist->FindBin(range.Stop());

        // try to autodedect maximum within fit range
        double maxx = range.Center();
        double maxy = -inf;
        for(int i=startbin; i<=stopbin; ++i) {
            const auto v = hist->GetBinContent(i);
            if(v > maxy) {
                maxy = v;
                maxx = hist->GetBinCenter(i);
            }
        }

        if(!isfinite(maxy))
           maxy = hist->GetMaximum();

        // linear background
        const math::LineFct bg({hist->GetBinCenter(startbin), hist->GetBinContent(startbin)},
                               {hist->GetBinCenter(stopbin),  hist->GetBinContent(stopbin)} );

        // amplitude
        func->SetParameter(0, maxy - bg(maxx));

        // x0
        func->SetParameter(1, maxx);

        // sigma
        func->SetParameter(2, 1.5);

        // pol1
        func->SetParameter(3, bg.b);
        func->SetParameter(4, bg.m);

        Sync();
    }

};

Energy::GUI_Banana::GUI_Banana(const string& basename,
                               OptionsPtr options,
                               Energy::CalibType& type,
                               const std::shared_ptr<DataManager>& calmgr,
                               const detector_ptr_t& detector,
                               const interval<double>& projectionrange,
                               const double proton_peak_mc_pos
                               ) :
    GUI_CalibType(basename, options, type, calmgr, detector),
    func(make_shared<FitProtonPeak>()),
    projection_range(projectionrange),
    proton_peak_mc(proton_peak_mc_pos),
    full_hist_name(
            options->Get<string>("HistogramPath", CalibModule_traits::GetName())
            + "/"
            + options->Get<string>("HistogramName", "Bananas"))
{

}

std::shared_ptr<TH1> Energy::GUI_Banana::GetHistogram(const WrapTFile& file) const
{
    return file.GetSharedHist<TH1>(full_hist_name);
}

void Energy::GUI_Banana::InitGUI(gui::ManagerWindow_traits* window)
{
    GUI_CalibType::InitGUI(window);
    window->AddNumberEntry("Chi2/NDF limit for autostop", AutoStopOnChi2);

    c_fit = window->AddCalCanvas();
    c_extra = window->AddCalCanvas();


    h_relative = new TH1D("h_relative","Relative change from previous gains",GetNumberOfChannels(),0,GetNumberOfChannels());
    h_relative->SetXTitle("Channel Number");
    h_relative->SetYTitle("Relative change / %");
}

gui::CalibModule_traits::DoFitReturn_t Energy::GUI_Banana::DoFit(const TH1& hist, unsigned ch)
{
    if(detector->IsIgnored(ch))
        return DoFitReturn_t::Skip;

    auto h_bananas = dynamic_cast<const TH3&>(hist);
    h_bananas.GetZaxis()->SetRange(ch+1,ch+1);
    banana = dynamic_cast<TH2D*>(h_bananas.Project3D("yx"));
    auto xaxis = banana->GetXaxis();
    h_projection = dynamic_cast<TH1D*>(banana->ProjectionY(
                                           "_py",
                                           xaxis->FindFixBin(projection_range.Start()),
                                           xaxis->FindFixBin(projection_range.Stop())
                                           )
                                       );

    if(h_projection->GetNbinsX() > 100) {
        const auto grp = int(std::ceil(h_projection->GetNbinsX()/100.0));
        h_projection->Rebin(grp);
    }

    // stop at empty histograms
    if(h_projection->GetEntries()==0)
        return DoFitReturn_t::Display;

    func->SetRange(interval<double>(0.5,6));
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

void Energy::GUI_Banana::DisplayFit()
{
    c_fit->Show(h_projection, func.get());

    c_extra->cd();
    banana->Draw("colz");
}

void Energy::GUI_Banana::StoreFit(unsigned channel)
{
    const double oldValue = previousValues[channel];

    const double protonpeak = func->GetPeakPosition();

    const double newValue = oldValue * proton_peak_mc / protonpeak;

    calibType.Values[channel] = newValue;

    const double relative_change = 100*(newValue/oldValue-1);

    LOG(INFO) << "Stored Ch=" << channel << ": ProtonPeak " << protonpeak
              << " MeV,  gain changed " << oldValue << " -> " << newValue
              << " (" << relative_change << " %)";


    // don't forget the fit parameters
    fitParameters[channel] = func->Save();

    h_relative->SetBinContent(channel+1, relative_change);

}

bool Energy::GUI_Banana::FinishSlice()
{
    c_extra->Clear();
    c_fit->Clear();

    c_fit->cd();
    h_relative->SetStats(false);
    h_relative->Draw("P");

    return true;
}


Energy::GUI_MIP::GUI_MIP(const string& basename,
                         OptionsPtr options,
                         Energy::CalibType& type,
                         const std::shared_ptr<DataManager>& calmgr,
                         const detector_ptr_t& detector,
                         const double peak_mc_pos
                         ) :
    GUI_CalibType(basename, options, type, calmgr, detector),
    func(make_shared<gui::FitLandauExpo>()),
    peak_mc(peak_mc_pos),
    full_hist_name(
            options->Get<string>("HistogramPath", CalibModule_traits::GetName())
            + "/"
            + options->Get<string>("HistogramName", "MIP"))
{

}

std::shared_ptr<TH1> Energy::GUI_MIP::GetHistogram(const WrapTFile& file) const
{
    return file.GetSharedHist<TH1>(full_hist_name);
}

void Energy::GUI_MIP::InitGUI(gui::ManagerWindow_traits* window)
{
    GUI_CalibType::InitGUI(window);
    window->AddNumberEntry("Chi2/NDF limit for autostop", AutoStopOnChi2);

    canvas = window->AddCalCanvas();

    h_peaks = new TH1D("h_peaks","Peak positions",GetNumberOfChannels(),0,GetNumberOfChannels());
    h_peaks->SetXTitle("Channel Number");
    h_peaks->SetYTitle("Minimum Ionizing Peak / MeV");
    h_relative = new TH1D("h_relative","Relative change from previous gains",GetNumberOfChannels(),0,GetNumberOfChannels());
    h_relative->SetXTitle("Channel Number");
    h_relative->SetYTitle("Relative change / %");
}

gui::CalibModule_traits::DoFitReturn_t Energy::GUI_MIP::DoFit(const TH1& hist, unsigned ch)
{
    if(detector->IsIgnored(ch))
        return DoFitReturn_t::Skip;

    auto hist2 = dynamic_cast<const TH2&>(hist);
    h_projection = hist2.ProjectionX("h_projection",ch+1,ch+1);

    // stop at empty histograms
    if(h_projection->GetEntries()==0)
        return DoFitReturn_t::Display;

    auto range = interval<double>(0.5,7);

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

    // try with defaults and signal fit
    func->SetDefaults(h_projection);
    func->SetRange(range);
    func->FitSignal(h_projection);

    if(fit_loop(5))
        return DoFitReturn_t::Next;


    // reached maximum retries without good chi2
    LOG(INFO) << "Chi2/dof = " << func->Chi2NDF();
    return DoFitReturn_t::Display;
}

void Energy::GUI_MIP::DisplayFit()
{
    canvas->Show(h_projection, func.get());
}

void Energy::GUI_MIP::StoreFit(unsigned channel)
{
    const double oldValue = previousValues[channel];
    const double peak = func->GetPeakPosition();
    const double newValue = oldValue * peak_mc / peak;

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

bool Energy::GUI_MIP::FinishSlice()
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


Energy::GUI_HEP::GUI_HEP(const string& basename,
                         OptionsPtr options,
                         Energy::CalibType& type,
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

std::shared_ptr<TH1> Energy::GUI_HEP::GetHistogram(const WrapTFile& file) const
{
    return file.GetSharedHist<TH1>(full_hist_name);
}

void Energy::GUI_HEP::InitGUI(gui::ManagerWindow_traits* window)
{
    GUI_CalibType::InitGUI(window);
    window->AddNumberEntry("Chi2/NDF limit for autostop", AutoStopOnChi2);

    canvas = window->AddCalCanvas();

    h_peaks = new TH1D("h_peaks","Peak positions",GetNumberOfChannels(),0,GetNumberOfChannels());
    h_peaks->SetXTitle("Channel Number");
    h_peaks->SetYTitle("High Energy Proton Peak / MeV");
    h_relative = new TH1D("h_relative","Relative change from previous gains",GetNumberOfChannels(),0,GetNumberOfChannels());
    h_relative->SetXTitle("Channel Number");
    h_relative->SetYTitle("Relative change / %");
}

gui::CalibModule_traits::DoFitReturn_t Energy::GUI_HEP::DoFit(const TH1& hist, unsigned ch)
{
    if(detector->IsIgnored(ch))
        return DoFitReturn_t::Skip;

    auto hist2 = dynamic_cast<const TH2&>(hist);
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

void Energy::GUI_HEP::DisplayFit()
{
    canvas->Show(h_projection, func.get());
}

void Energy::GUI_HEP::StoreFit(unsigned channel)
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

bool Energy::GUI_HEP::FinishSlice()
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
