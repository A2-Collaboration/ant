#include "Time.h"
#include "analysis/plot/HistogramFactories.h"
#include "analysis/data/Event.h"
#include "analysis/utils/combinatorics.h"

#include "calibration/gui/CalCanvas.h"
#include "calibration/fitfunctions/FitFunction.h" // PeakingFitFunction

#include "DataManager.h"

#include "tree/TDetectorRead.h"
#include "tree/TCalibrationData.h"

#include "base/Logger.h"

#include "TH1D.h"
#include "TH2D.h"

#include <cstdint>

using namespace std;
using namespace ant;
using namespace ant::calibration;
using namespace ant::analysis;
using namespace ant::analysis::data;

Time::Time(const std::shared_ptr<Detector_t>& detector, const std::shared_ptr<DataManager>& CalibrationManager,
           Calibration::Converter::ptr_t converter,
           double defaultOffset,
           shared_ptr<gui::PeakingFitFunktion> FitFunction,
           const interval<double>& timeWindow, // default {-inf, inf}
           const double defaultGain, // default gain is 1.0
           const std::vector< TKeyValue<double> >& gains) :
    Calibration::Module(
        std_ext::formatter()
        << Detector_t::ToString(detector->Type)
        << "_Time"
           ),
    Detector(detector),
    calibrationManager(CalibrationManager),
    Converter(move(converter)),
    TimeWindow(timeWindow),
    fitFunction(FitFunction),
    DefaultOffset(defaultOffset),
    Offsets(),
    DefaultGain(defaultGain),
    Gains()
{
    if(Converter==nullptr)
        throw std::runtime_error("Given converter should not be nullptr");

    // fill a gain vector from given key-value pairs
    // for faster access (if some are given at all)
    if(gains.empty())
        return;
    unsigned maxkey = 0;
    for(const auto& gain : gains)
        maxkey = gain.Key>maxkey ? gain.Key : maxkey;
    Gains.resize(maxkey+1, DefaultGain);
    for(const auto& gain : gains)
        Gains[gain.Key] = gain.Value;
}

std::vector<std::list<TID> > Time::GetChangePoints() const {
    vector<list<TID>> changePointLists;
    changePointLists.emplace_back(calibrationManager->GetChangePoints(GetName()));
    return changePointLists;
}

void Time::Update(size_t, const TID& id)
{
    TCalibrationData cdata;
    if(calibrationManager->GetData(GetName(), id, cdata))
    {
        for (const auto& val: cdata.Data) {
            if(Offsets.size()<val.Key+1)
                Offsets.resize(val.Key+1);
            Offsets[val.Key] = val.Value;
        }
    }
    else {
        LOG(WARNING) << "No calibration data found for offsets"
                     << " at changepoint TID=" << id << ", using previous values";
    }
}

void Time::ApplyTo(const readhits_t& hits, extrahits_t&)
{
    const auto& dethits = hits.get_item(Detector->Type);

    // now calibrate the Times (ignore any other kind of hits)
    for(TDetectorReadHit* dethit : dethits) {
        if(dethit->GetChannelType() != Channel_t::Type_t::Timing)
            continue;

        // the Converter is smart enough to account for reference Times!
        const auto& values = Converter->Convert(dethit->RawData);
        dethit->Values.reserve(values.size());

        // apply gain/offset to each of the values (might be multihit)
        for(double value : values) {
            if(Gains.empty())
                value *= DefaultGain;
            else
                value *= Gains[dethit->Channel];

            if(Offsets.empty())
                value -= DefaultOffset;
            else
                value -= Offsets[dethit->Channel];

            if(!TimeWindow.Contains(value))
            {
                VLOG(9) << "Discarding hit in channel " << dethit->Channel << ", which is outside time window.";
                continue;
            }

            dethit->Values.push_back(value);
        }
    }
}

Time::ThePhysics::ThePhysics(const string& name, const string& histName,
                             const std::shared_ptr<Detector_t>& theDetector):
    Physics(name),
    detector(theDetector)
{
    string detectorName(Detector_t::ToString(detector->Type));
    hTime = HistFac.makeTH2D( detectorName + string(" - Time"),
                              "time [ns]",
                              detectorName + "-channel",
                              BinSettings(1000,-400,400),
                              BinSettings(detector->GetNChannels()),
                              histName
                              );
}


void Time::ThePhysics::ProcessEvent(const Event& event)
{
    //handle Tagger differently
    if (detector->Type == Detector_t::Type_t::EPT)
    {
        for (const auto& tHit: event.Reconstructed().TaggerHits())
            hTime->Fill(tHit->Time(),tHit->Channel());
        return;
    }

    for ( const auto& cand: event.Reconstructed().Candidates())
        for (const auto& cluster: cand->Clusters)
            if (cluster.Detector == detector->Type)
                hTime->Fill(cluster.Time,cluster.CentralElement);

}

void Time::ThePhysics::Finish()
{
}

void Time::ThePhysics::ShowResult()
{
    canvas(GetName())  << drawoption("colz") << hTime << endc;
}


Time::TheGUI::TheGUI(const string& name,
                     const std::shared_ptr<Detector_t>& theDetector,
                     const std::shared_ptr<DataManager>& cDataManager,
                     double DefaultOffset,
                     const std::vector<double>& Offsets,
                     const shared_ptr<gui::PeakingFitFunktion> FitFunction):
    gui::Manager_traits(name),
    detector(theDetector),
    calmgr(cDataManager),
    defaultOffset(DefaultOffset),
    offsets(Offsets),
    theCanvas(nullptr),
    times(nullptr),
    fitFunction(FitFunction)
{
}

void Time::TheGUI::InitGUI(gui::ManagerWindow_traits* window)
{
    theCanvas = window->AddCalCanvas();
    times = new TH1D("times","Times",
                     1000, -400 ,400);
    times->SetXTitle("time [ns]");
    times->SetYTitle("#");
    timePeaks = new TH1D("timePeaks","Time Peaks",
                         GetNumberOfChannels(), 0, GetNumberOfChannels());
    timePeaks->SetXTitle("Channel");
    timePeaks->SetYTitle("Peak position [ns]");
}

void Time::TheGUI::StartRange(const interval<TID>& range)
{
    offsets.resize(GetNumberOfChannels(),defaultOffset);
    TCalibrationData cdata;
    if (calmgr->GetData(GetName(),range.Start(),cdata))
    {
        for (const TKeyValue<double>& entry: cdata.Data)
            offsets[entry.Key] = entry.Value;
        for (const TKeyValue<vector<double>>& entry: cdata.FitParameters)
            fitParams.insert(make_pair(entry.Key,entry.Value));
        LOG(INFO) << GetName() << ": Loaded previous values from database";
    }
    else
        LOG(INFO) << GetName() << ": No previous offsets found, built from default offset";

    previousOffsets = offsets;
}

gui::Manager_traits::DoFitReturn_t Time::TheGUI::DoFit(TH1* hist, unsigned channel)
{
    if (detector->IsIgnored(channel))
        return gui::Manager_traits::DoFitReturn_t::Skip;


    TH2* hist2 = dynamic_cast<TH2*>(hist);

    times = hist2->ProjectionX("",channel+1,channel+1);

    fitFunction->SetDefaults(times);
    const auto it_fit_param = fitParams.find(channel);
    if(it_fit_param != fitParams.end()) {
        VLOG(5) << "Loading previous fit parameters for channel " << channel;
        fitFunction->Load(it_fit_param->second);
    }

    fitFunction->Fit(times);

    /// \todo auto-stop if fit fails

    // do not show something, goto next channel
    return DoFitReturn_t::Next;

}

void Time::TheGUI::DisplayFit()
{
    theCanvas->Show(times, fitFunction.get());
}

void Time::TheGUI::StoreFit(unsigned channel)
{
    const double oldOffset = previousOffsets[channel];
    const double timePeak = fitFunction->GetPeakPosition();

    timePeaks->Fill(channel,timePeak);

    // the timePeak should be zero, so this gives directly
    // the value to change the offset
    const double newOffset = oldOffset - timePeak;

    offsets[channel] = newOffset;

    const double relative_change = 100*(newOffset/oldOffset-1);

    LOG(INFO) << "Stored Ch=" << channel << ": PeakPosition " << timePeak
              << " ns,  offset changed " << oldOffset << " -> " << newOffset
              << " (" << relative_change << " %)";


    // don't forget the fit parameters
    fitParams[channel] = fitFunction->Save();

    theCanvas->Clear();
    theCanvas->Update();
}

bool Time::TheGUI::FinishRange()
{
    theCanvas->Clear();
    theCanvas->cd();
    timePeaks->SetStats(false);
    timePeaks->SetMarkerStyle(20);

    timePeaks->Draw("P");

    theCanvas->Update();

    return true;
}

void Time::TheGUI::StoreFinishRange(const interval<TID>& range)
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

    LOG(INFO) << "Added " << cdata;
}
