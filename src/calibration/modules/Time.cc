#include "Time.h"


#include "calibration/gui/CalCanvas.h"
#include "calibration/fitfunctions/FitFunction.h" // PeakingFitFunction

#include "DataManager.h"

#include "tree/TDetectorReadHit.h"
#include "tree/TCalibrationData.h"

#include "base/Logger.h"
#include "base/PlotExt.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"

#include <cstdint>

using namespace std;
using namespace ant;
using namespace ant::calibration;

Time::Time(const std::shared_ptr<Detector_t>& detector, const std::shared_ptr<DataManager>& CalibrationManager,
           Calibration::Converter::ptr_t converter,
           double defaultOffset,
           shared_ptr<gui::PeakingFitFunction> FitFunction,
           const interval<double>& timeWindow, // default {-inf, inf}
           const double defaultGain) :
    Calibration::Module(
        std_ext::formatter()
        << Detector_t::ToString(detector->Type)
        << "_Time"
           ),
    Detector(detector),
    calibrationManager(CalibrationManager),
    Converters(detector->GetNChannels(), converter),
    TimeWindows(detector->GetNChannels(), timeWindow),
    fitFunction(FitFunction),
    DefaultOffsets(detector->GetNChannels(), defaultOffset),
    Offsets(),
    DefaultGains(detector->GetNChannels(), defaultGain),
    Gains()
{
}

std::list<Updateable_traits::Loader_t> Time::GetLoaders()
{
    return {
      [this] (const TID& currPoint, TID& nextChangePoint) {
            TCalibrationData cdata;
            if(calibrationManager->GetData(GetName(), currPoint, cdata, nextChangePoint))
            {
                for (const auto& val: cdata.Data) {
                    if(Offsets.size()<val.Key+1)
                        Offsets.resize(val.Key+1);
                    Offsets[val.Key] = val.Value;
                }
            }
            else {
                LOG_IF(!Offsets.empty(), WARNING) << "No calibration data found for offsets"
                                                  << " at changepoint TID="
                                                  << currPoint << ", using default values";
                Offsets.resize(0);
            }
        }
    };
}

void Time::UpdatedTIDFlags(const TID& id)
{
    IsMC = id.isSet(TID::Flags_t::MC);
}

void Time::GetGUIs(std::list<std::unique_ptr<gui::CalibModule_traits> >& guis, OptionsPtr) {
    guis.emplace_back(std_ext::make_unique<TheGUI>(
                          GetName(),
                          Detector,
                          calibrationManager,
                          DefaultOffsets,
                          Offsets,
                          fitFunction
                          ));
}

void Time::ApplyTo(const readhits_t& hits)
{
    /// \bug MC could also potentially be calibrated
    if(IsMC)
        return;

    auto& dethits = hits.get_item(Detector->Type);

    // now calibrate the Times (ignore any other kind of hits)
    for(TDetectorReadHit& dethit : dethits) {
        if(dethit.ChannelType != Channel_t::Type_t::Timing)
            continue;

        // clear possible previous reads
        dethit.Values.resize(0);

        // the Converter is smart enough to account for reference times
        // by (possibly) being itself a reconstruction hook and searching for it
        const auto& converted = Converters[dethit.Channel]->Convert(dethit.RawData);

        // apply gain/offset to each of the values (might be multihit)
        for(const double& conv : converted) {
            TDetectorReadHit::Value_t value(conv);

            if(Gains.empty())
                value.Calibrated *= DefaultGains[dethit.Channel];
            else
                value.Calibrated *= Gains[dethit.Channel];

            if(Offsets.empty())
                value.Calibrated -= DefaultOffsets[dethit.Channel];
            else
                value.Calibrated -= Offsets[dethit.Channel];

            if(!TimeWindows[dethit.Channel].Contains(value.Calibrated))
            {
                VLOG(9) << "Discarding hit in channel " << dethit.Channel << ", which is outside time window.";
                continue;
            }

            dethit.Values.emplace_back(move(value));
        }
    }
}

Time::TheGUI::TheGUI(const string& name,
                     const std::shared_ptr<Detector_t>& theDetector,
                     const std::shared_ptr<DataManager>& cDataManager,
                     const std::vector<double>& DefaultOffsets,
                     const std::vector<double>& Offsets,
                     const shared_ptr<gui::PeakingFitFunction> FitFunction):
    gui::CalibModule_traits(name),
    detector(theDetector),
    calmgr(cDataManager),
    defaultOffsets(DefaultOffsets),
    offsets(Offsets),
    theCanvas(nullptr),
    times(nullptr),
    fitFunction(FitFunction)
{
    // especially for histogram with low statistics,
    // this improves the fit result (does not fit weird noise spikes anymore)
    fitFunction->SetAdditionalFitArgs("W");
}

shared_ptr<TH1> Time::TheGUI::GetHistogram(const WrapTFile& file) const {
    return file.GetSharedHist<TH1>(GetName()+"/Time");
}

unsigned Time::TheGUI::GetNumberOfChannels() const
{
    return detector->GetNChannels();
}

void Time::TheGUI::InitGUI(gui::ManagerWindow_traits& window)
{
    window.AddCheckBox("Ignore prev fit params", IgnorePreviousFitParameters);

    window.AddNumberEntry("Chi2/NDF limit for autostop", AutoStopOnChi2);
    window.AddNumberEntry("Max Peakposition difference", AutoStopOnPeakPos);
    window.AddNumberEntry("Stop at Channel", AutoStopAtChannel);
    window.AddNumberEntry("time in [ - t_0 , t_0 ]", HardTimeCut);
    window.AddNumberEntry("Rebin", Rebin);

    window.AddCheckBox("Skip empty channels", SkipEmptyChannels);

    theCanvas = window.AddCalCanvas();
    times = new TH1D("times","Times",
                     1000, -400 ,400);
    times->SetXTitle("time [ns]");
    times->SetYTitle("#");
    timePeaks = new TH1D("timePeaks","Time Peaks",
                         GetNumberOfChannels(), 0, GetNumberOfChannels());
    timePeaks->SetXTitle("Channel");
    timePeaks->SetYTitle("Peak position [ns]");
}

void Time::TheGUI::StartSlice(const interval<TID>& range)
{
    offsets = defaultOffsets;
    TCalibrationData cdata;
    if (calmgr->GetData(GetName(),range.Start(),cdata))
    {
        for (const TKeyValue<double>& entry: cdata.Data)
            offsets[entry.Key] = entry.Value;
        for (const TKeyValue<vector<double>>& entry: cdata.FitParameters)
            fitParams[entry.Key] = entry.Value;
        LOG(INFO) << GetName() << ": Loaded previous values from database";
    }
    else
        LOG(INFO) << GetName() << ": No previous offsets found, built from default offset";

    // remember the previous offsets
    previousOffsets = offsets;
}

gui::CalibModule_traits::DoFitReturn_t Time::TheGUI::DoFit(const TH1& hist, unsigned channel)
{
    if (detector->IsIgnored(channel))
        return gui::CalibModule_traits::DoFitReturn_t::Skip;


    auto& hist2 = dynamic_cast<const TH2&>(hist);

    times = hist2.ProjectionX("times",channel+1,channel+1);
    if(Rebin>1.0)
        times = times->Rebin(int(Rebin));


    if(times->GetEntries() == 0 && SkipEmptyChannels) {
        channelWasEmpty = true;
        return DoFitReturn_t::Next;
    } else {
        channelWasEmpty = false;
    }

    if (HardTimeCut > 0 )
        times->GetXaxis()->SetRangeUser(-fabs(HardTimeCut),fabs(HardTimeCut));

    fitFunction->SetDefaults(times);
    const auto it_fit_param = fitParams.find(channel);
    if(it_fit_param != fitParams.end() && !IgnorePreviousFitParameters) {
        VLOG(5) << "Loading previous fit parameters for channel " << channel;
        fitFunction->Load(it_fit_param->second);
    }

    const auto maximum = GetMaxPos(times);
    bool chi2OK = false;
    bool PeakPosOK = false;
    size_t retries = 5;
    do {
        fitFunction->Fit(times);
        VLOG(5) << "Chi2/dof = " << fitFunction->Chi2NDF();

        chi2OK = fitFunction->Chi2NDF() < AutoStopOnChi2 ;
        PeakPosOK = fabs(maximum - fitFunction->GetPeakPosition()) < AutoStopOnPeakPos ;
        if( chi2OK && PeakPosOK )  {
            break;
        }
        retries--;
    }
    while(retries>0);

    if(int(AutoStopAtChannel) == int(channel))
        return DoFitReturn_t::Display;

    if( chi2OK && PeakPosOK )
        return DoFitReturn_t::Next;

    // reached maximum retries without good chi2

    LOG(INFO) << "Stopped automode" ;
    if (!chi2OK) LOG(INFO)    << " -> Chi2/dof = " << fitFunction->Chi2NDF();

    if (!PeakPosOK) LOG(INFO) << " -> Distance Max to PeakPos : " << maximum << " - " <<  fitFunction->GetPeakPosition()
                              << " = " << fabs(maximum - fitFunction->GetPeakPosition());

    return DoFitReturn_t::Display;

}

void Time::TheGUI::DisplayFit()
{
    theCanvas->Show(times, fitFunction.get());
}

void Time::TheGUI::StoreFit(unsigned channel)
{
    const double oldOffset = previousOffsets[channel];
    const double timePeak = !channelWasEmpty ? fitFunction->GetPeakPosition() : 0.0 ;

    timePeaks->SetBinContent(channel+1,timePeak);

    // the timePeak should be zero, so this gives directly
    // the value to change the offset
    const double newOffset = oldOffset + timePeak;

    offsets[channel] = newOffset;

    const double relative_change = 100*(newOffset/oldOffset-1);

    LOG(INFO) << "Stored Ch=" << channel << ": PeakPosition " << timePeak
              << " ns,  offset changed " << oldOffset << " -> " << newOffset
              << " (" << relative_change << " %)";


    if(!channelWasEmpty) {
        // don't forget the fit parameters
        fitParams[channel] = fitFunction->Save();
    }

    theCanvas->Clear();
    theCanvas->Update();
}

bool Time::TheGUI::FinishSlice()
{
    theCanvas->Clear();
    theCanvas->cd();
    timePeaks->SetStats(false);
    timePeaks->SetMarkerStyle(20);

    timePeaks->Draw("P");

    theCanvas->Update();

    return true;
}

void Time::TheGUI::StoreFinishSlice(const interval<TID>& range)
{
    theCanvas->Clear();
    theCanvas->Update();

    TCalibrationData cdata(
                GetName(),
                range.Start(),
                range.Stop()
                );


    // fill data
    for(unsigned ch=0;ch<offsets.size();ch++) {
        cdata.Data.emplace_back(ch, offsets[ch]);
    }

    // fill fit parameters (if any)
    for(const auto& it_map : fitParams) {
        const unsigned ch = it_map.first;
        const vector<double>& params = it_map.second;
        cdata.FitParameters.emplace_back(ch, params);
    }

    calmgr->Add(cdata);
}

void gui::CBPeakFunction::SetDefaults(TH1* hist)
{
    if(hist) {
        func->SetParameter(0,hist->GetMaximum());
        const double max_pos = hist->GetXaxis()->GetBinCenter(hist->GetMaximumBin());
        func->SetParameter(1,max_pos);
        const double sigma = 2;
        func->SetParameter(2, sigma);
        SetRange({max_pos-40*sigma, max_pos+40*sigma});
    } else {
        SetRange({0,200});
        func->SetParameter(0,0.8);
        func->SetParameter(1,100);
        func->SetParameter(2,20);
    }

}
