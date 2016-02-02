#include "Energy.h"

#include "calibration/DataManager.h"
#include "calibration/gui/CalCanvas.h"
#include "calibration/fitfunctions/FitGausPol1.h"

#include "tree/TCalibrationData.h"
#include "tree/TDetectorReadHit.h"

#include "base/Logger.h"

#include "TH2.h"
#include "TH3.h"
#include "TF1.h"

#include <cstdint>
#include <ctime>
#include <algorithm>
#include <sstream>
#include <vector>
#include <list>

using namespace std;
using namespace ant;
using namespace ant::calibration;

Energy::Energy(Detector_t::Type_t detectorType,
               std::shared_ptr<DataManager> calmgr,
               Calibration::Converter::ptr_t converter,
               double defaultPedestal,
               double defaultGain,
               double defaultThreshold,
               double defaultRelativeGain,
               Channel_t::Type_t channelType) :
    Calibration::Module(
        std_ext::formatter()
        << Detector_t::ToString(detectorType)
        << "_"
        << ( channelType == Channel_t::Type_t::IntegralShort ? "Short" : "" )
        << "Energy"
           ),
    DetectorType(detectorType),
    ChannelType(channelType),
    calibrationManager(calmgr),
    Converter(move(converter)),
    Pedestals(defaultPedestal, "Pedestals"),
    Gains(defaultGain, "Gains", "ggIM"),
    Thresholds(defaultThreshold, "Thresholds"),
    RelativeGains(defaultRelativeGain, "RelativeGains", "ggIM")
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
    for(TDetectorReadHit* dethit : dethits) {
        if(dethit->ChannelType != ChannelType)
            continue;



        // Values might already be filled
        // (for example by previous calibration run, or A2Geant unpacker),
        // then we apply the threshold and the relative gain only
        std::vector<double> all_values;

        // prefer RawData if available
        if(!dethit->RawData.empty()) {
            // convert to not-so-raw values (still not MeV scale though)
            dethit->Converted = Converter->Convert(dethit->RawData);

            // apply pedestal/gain to each of the values (might be multihit)
            for(double value : dethit->Converted) {
                value -= Pedestals.Get(dethit->Channel);
                value *= Gains.Get(dethit->Channel);
                all_values.push_back(value);
            }

        }
        else {
            // maybe the values are already filled
            all_values = dethit->Values;
        }

        // always apply the threshold cut and the relative gains
        dethit->Values.resize(0);
        dethit->Values.reserve(all_values.size());

        for(double value : all_values) {
            value *= RelativeGains.Get(dethit->Channel);

            const double threshold = Thresholds.Get(dethit->Channel);
            if(value<threshold)
                continue;

            // only add if it passes the threshold
            dethit->Values.push_back(value);
        }

    }
}

double Energy::CalibType::Get(unsigned channel) const {
    if(Values.empty()) {
        if(DefaultValues.empty()) {
            return DefaultValue;
        }
        else {
            return DefaultValues[channel];
        }
    }
    else {
        return Values[channel];
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

Energy::GUI_CalibType::GUI_CalibType(const string& basename, CalibType& type,
                                     const shared_ptr<DataManager>& calmgr,
                                     const shared_ptr<Detector_t>& detector_,
                                     Calibration::AddMode_t mode) :
    gui::CalibModule_traits(basename),
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
    return file.GetSharedHist<TH1>(CalibModule_traits::GetName()+"/"+calibType.HistogramName);
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
    std::vector<double>& values = calibType.Values;
    values.resize(GetNumberOfChannels(), calibType.DefaultValue);

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

Energy::GUI_Pedestals::GUI_Pedestals(
        const string& basename,
        CalibType& type,
        const std::shared_ptr<DataManager>& calmgr,
        const std::shared_ptr<Detector_t>& detector,
        shared_ptr<gui::PeakingFitFunction> fitfunction) :
    GUI_CalibType(basename, type, calmgr, detector, Calibration::AddMode_t::RightOpen),
    func(fitfunction)
{

}

void Energy::GUI_Pedestals::InitGUI(gui::ManagerWindow_traits* window)
{
    GUI_CalibType::InitGUI(window);
    canvas = window->AddCalCanvas();
}

gui::CalibModule_traits::DoFitReturn_t Energy::GUI_Pedestals::DoFit(TH1* hist, unsigned channel)
{
    if(detector->IsIgnored(channel))
        return DoFitReturn_t::Skip;

    TH2* hist2 = dynamic_cast<TH2*>(hist);

    h_projection = hist2->ProjectionX("h_projection",channel+1,channel+1);

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

        // amplitude
        func->SetParameter(0, 0.5*hist->GetMaximum());

        // x0
        func->SetParameter(1, 8);

        // sigma
        func->SetParameter(2, 1.5);

        // pol1
        func->SetParameter(3, 0.1*hist->GetMaximum());
        func->SetParameter(4, 0);

        Sync();
    }

};

Energy::GUI_Banana::GUI_Banana(const string& basename,
                                 Energy::CalibType& type,
                                 const std::shared_ptr<DataManager>& calmgr,
                                 const std::shared_ptr<Detector_t>& detector) :
    GUI_CalibType(basename, type, calmgr, detector),
    func(make_shared<FitProtonPeak>())
{

}

std::shared_ptr<TH1> Energy::GUI_Banana::GetHistogram(const WrapTFile& file) const
{
    return file.GetSharedHist<TH1>(CalibModule_traits::GetName()+"/Bananas");
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

gui::CalibModule_traits::DoFitReturn_t Energy::GUI_Banana::DoFit(TH1* hist, unsigned ch)
{
    if(detector->IsIgnored(ch))
        return DoFitReturn_t::Skip;

    TH3* h_bananas = dynamic_cast<TH3*>(hist);
    h_bananas->GetZaxis()->SetRange(ch+1,ch+1);
    banana = dynamic_cast<TH2D*>(h_bananas->Project3D("yx"));
    auto xaxis = banana->GetXaxis();
    h_projection = dynamic_cast<TH1D*>(banana->ProjectionY(
                                           "_py",
                                           xaxis->FindFixBin(130),
                                           xaxis->FindFixBin(150)
                                           )
                                       );

    // stop at empty histograms
    if(h_projection->GetEntries()==0)
        return DoFitReturn_t::Display;

    func->SetRange(interval<double>(5,13));
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
    /// \todo obtain convergenceFactor and pi0mass from config or database
    const double protonMC = 8.0;
    const double protonpeak = func->GetPeakPosition();

    const double newValue = oldValue * protonMC/protonpeak;

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
