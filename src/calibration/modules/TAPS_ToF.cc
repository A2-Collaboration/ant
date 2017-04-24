#include "TAPS_ToF.h"

#include "calibration/gui/CalCanvas.h"
#include "fitfunctions/FitLandauPol0.h"

#include "DataManager.h"

#include "tree/TCalibrationData.h"
#include "expconfig/detectors/TAPS.h"

#include "base/Logger.h"

#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"

#include <cstdint>

using namespace std;
using namespace ant;
using namespace ant::calibration;

TAPS_ToF::TAPS_ToF(const std::shared_ptr<expconfig::detector::TAPS>& detector,
               const std::shared_ptr<DataManager>& CalibrationManager) :
    Calibration::Module("TAPS_ToF_Offsets"),
    Detector(detector),
    calibrationManager(CalibrationManager)
{

}

std::list<Updateable_traits::Loader_t> TAPS_ToF::GetLoaders()
{
    return {
      [this] (const TID& currPoint, TID& nextChangePoint) {
            TCalibrationData cdata;
            if(calibrationManager->GetData(GetName(), currPoint, cdata, nextChangePoint))
            {
                for (const auto& val: cdata.Data) {
                    Detector->SetToFOffset(val.Key, val.Value);
                }
            }
        }
    };
}



void TAPS_ToF::GetGUIs(std::list<std::unique_ptr<gui::CalibModule_traits> >& guis, ant::OptionsPtr) {
    guis.emplace_back(std_ext::make_unique<TheGUI>(
                          GetName(),
                          Detector,
                          calibrationManager
                          ));
}

std::vector<string> TAPS_ToF::GetPhysicsModules() const {
    return {"TAPS_Time"}; // use Time_ToF histogram from there
}

TAPS_ToF::TheGUI::TheGUI(const string& name,
                         const std::shared_ptr<Detector_t>& detector_,
                         const std::shared_ptr<DataManager>& cDataManager) :
    gui::CalibModule_traits(name),
    detector(detector_),
    calmgr(cDataManager),
    fitFunction(make_shared<gui::FitLandauPol0>())
{
}

shared_ptr<TH1> TAPS_ToF::TheGUI::GetHistogram(const WrapTFile& file) const {
    return file.GetSharedHist<TH1>("TAPS_Time/hTimeToTriggerRef");
}

unsigned TAPS_ToF::TheGUI::GetNumberOfChannels() const
{
    return detector->GetNChannels();
}

void TAPS_ToF::TheGUI::InitGUI(gui::ManagerWindow_traits& window)
{
    window.AddCheckBox("Ignore prev fit params", IgnorePreviousFitParameters);

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

void TAPS_ToF::TheGUI::StartSlice(const interval<TID>& range)
{
    offsets.resize(0);
    offsets.resize(GetNumberOfChannels(), 0);
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
        LOG(INFO) << GetName() << ": No previous offsets found, built from default offset 0";

    // remember the previous offsets
    previousOffsets = offsets;
}

gui::CalibModule_traits::DoFitReturn_t TAPS_ToF::TheGUI::DoFit(const TH1& hist, unsigned channel)
{
    if (detector->IsIgnored(channel))
        return gui::CalibModule_traits::DoFitReturn_t::Skip;


    auto& hist2 = dynamic_cast<const TH2&>(hist);

    times = hist2.ProjectionX("times",channel+1,channel+1);

    fitFunction->SetDefaults(times);
    const auto it_fit_param = fitParams.find(channel);
    if(it_fit_param != fitParams.end() && !IgnorePreviousFitParameters) {
        VLOG(5) << "Loading previous fit parameters for channel " << channel;
        fitFunction->Load(it_fit_param->second);
    }

    size_t retries = 5;
    do {
        fitFunction->Fit(times);
        VLOG(5) << "Chi2/dof = " << fitFunction->Chi2NDF();
        if(fitFunction->Chi2NDF() < 3000.0) {
            return DoFitReturn_t::Next;
        }
        retries--;
    }
    while(retries>0);

    // reached maximum retries without good chi2
    LOG(INFO) << "Chi2/dof = " << fitFunction->Chi2NDF();
    return DoFitReturn_t::Display;

}

void TAPS_ToF::TheGUI::DisplayFit()
{
    theCanvas->Show(times, fitFunction.get());
}

void TAPS_ToF::TheGUI::StoreFit(unsigned channel)
{
    const double oldOffset = previousOffsets[channel];
    const double timePeak = fitFunction->GetPeakPosition();

    timePeaks->SetBinContent(channel+1,timePeak);

    // the timePeak should be zero, so this gives directly
    // the value to change the offset
    const double newOffset = oldOffset + timePeak;

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

bool TAPS_ToF::TheGUI::FinishSlice()
{
    theCanvas->Clear();
    theCanvas->cd();
    timePeaks->SetStats(false);
    timePeaks->SetMarkerStyle(20);

    timePeaks->Draw("P");

    theCanvas->Update();

    return true;
}

void TAPS_ToF::TheGUI::StoreFinishSlice(const interval<TID>& range)
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
