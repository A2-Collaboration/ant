#include "CB_TimeWalk.h"

#include "DataManager.h"
#include "gui/CalCanvas.h"
#include "fitfunctions/FitTimewalk.h"
#include "tree/TCalibrationData.h"

#include "expconfig/detectors/CB.h"
#include "base/Logger.h"

#include "TGraph.h"
#include "TFitResult.h"
#include "TObjArray.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3.h"


#include <tree/TCluster.h>

#include <limits>
#include <cmath>

using namespace ant;
using namespace ant::calibration;
using namespace std;

CB_TimeWalk::CB_TimeWalk(const shared_ptr<expconfig::detector::CB>& cb,
        const shared_ptr<DataManager>& calmgr,
                         const interval<double>& timeWindow,
                         double badTDC_EnergyThreshold) :
    Module("CB_TimeWalk"),
    cb_detector(cb),
    calibrationManager(calmgr),
    TimeWindow(timeWindow),
    BadTDC_EnergyThreshold(badTDC_EnergyThreshold)
{
    for(unsigned ch=0;ch<cb_detector->GetNChannels();ch++) {
        timewalks.emplace_back(make_shared<gui::FitTimewalk>());
        timewalks.back()->SetDefaults(nullptr);
    }
}

CB_TimeWalk::~CB_TimeWalk()
{
}

void CB_TimeWalk::ApplyTo(clusterhits_t& sorted_clusterhits)
{
    if(IsMC)
        return;

    // search for CB clusters
    const auto it_sorted_clusterhits = sorted_clusterhits.find(Detector_t::Type_t::CB);
    if(it_sorted_clusterhits == sorted_clusterhits.end())
        return;

    auto& clusterhits = it_sorted_clusterhits->second;

    auto it_clusterhit = clusterhits.begin();

    while(it_clusterhit != clusterhits.end()) {
        TClusterHit& clusterhit = *it_clusterhit;

        // check if this hit belongs to a bad TDC chnanel
        if(cb_detector->HasElementFlags(clusterhit.Channel, Detector_t::ElementFlag_t::BadTDC)) {
            // if energy is above threhshold, keep the clusterhit
            // otherwise delete it
            if(clusterhit.Energy > BadTDC_EnergyThreshold)
                ++it_clusterhit;
            else
                it_clusterhit = clusterhits.erase(it_clusterhit);
            continue;
        }

        // do timewalk correction
        // use uncalibrated energy for that
        // to stay independent of energy calibration

        double raw_energy = std_ext::NaN;
        for(auto& datum : clusterhit.Data) {
            if(datum.Type == Channel_t::Type_t::Integral) {
                raw_energy = datum.Value.Uncalibrated;
                break;
            }
        }

        auto deltaT = timewalks[clusterhit.Channel]->Eval(raw_energy);
        if(isfinite(deltaT) && deltaT > -10 && deltaT < 80)
            clusterhit.Time -= deltaT;
        else if(isfinite(deltaT))
            LOG_N_TIMES(100, WARNING) << "(max 100 times) TimeWalk="
                                      << deltaT << " ignored for E=" << clusterhit.Energy
                                      << " ch=" << clusterhit.Channel;

        // get rid of clusterhit if outside timewindow
        if(std::isfinite(clusterhit.Time) && !TimeWindow.Contains(clusterhit.Time))
            it_clusterhit = clusterhits.erase(it_clusterhit);
        else
            ++it_clusterhit;
    }
}

void CB_TimeWalk::GetGUIs(list<unique_ptr<gui::CalibModule_traits> >& guis, OptionsPtr) {
    guis.emplace_back(std_ext::make_unique<TheGUI>(GetName(), calibrationManager, cb_detector, timewalks));
}


std::list<Updateable_traits::Loader_t> CB_TimeWalk::GetLoaders()
{
    return {
        [this] (const TID& currPoint, TID& nextChangePoint) {
            TCalibrationData cdata;
            if(!calibrationManager->GetData(GetName(), currPoint, cdata, nextChangePoint))
                return;
            for(const TKeyValue<vector<double>>& kv : cdata.FitParameters) {
                if(kv.Key>=timewalks.size()) {
                    LOG(ERROR) << "Ignoring too large key=" << kv.Key;
                    continue;
                }
                timewalks[kv.Key]->Load(kv.Value);
            }
        }
    };
}

void CB_TimeWalk::UpdatedTIDFlags(const TID& id)
{
    IsMC = id.isSet(TID::Flags_t::MC);
}

CB_TimeWalk::TheGUI::TheGUI(const string& basename,
                            const shared_ptr<DataManager>& calmgr,
                            const shared_ptr<expconfig::detector::CB>& cb,
                            std::vector< std::shared_ptr<gui::FitTimewalk> >& timewalks_) :
    gui::CalibModule_traits(basename),
    calibrationManager(calmgr),
    cb_detector(cb),
    timewalks(timewalks_)
{
}

shared_ptr<TH1> CB_TimeWalk::TheGUI::GetHistogram(const WrapTFile& file) const
{
    return file.GetSharedHist<TH1>(GetName()+"/timewalk");
}

unsigned CB_TimeWalk::TheGUI::GetNumberOfChannels() const
{
    return cb_detector->GetNChannels();
}

void CB_TimeWalk::TheGUI::InitGUI(gui::ManagerWindow_traits* window)
{
    c_fit = window->AddCalCanvas();
    c_extra = window->AddCalCanvas();
}

void CB_TimeWalk::TheGUI::StartSlice(const interval<TID>& range)
{

    TCalibrationData cdata;
    if(!calibrationManager->GetData(GetName(), range.Start(), cdata)) {
        LOG(INFO) << " No previous data found";
        return;
    }
    for(const TKeyValue<vector<double>>& kv : cdata.FitParameters) {
        if(kv.Key>=timewalks.size()) {
            LOG(ERROR) << "Ignoring too large key=" << kv.Key;
            continue;
        }
        timewalks[kv.Key]->Load(kv.Value);
    }
}



gui::CalibModule_traits::DoFitReturn_t CB_TimeWalk::TheGUI::DoFit(TH1* hist, unsigned ch)
{
    if(cb_detector->IsIgnored(ch))
        return DoFitReturn_t::Skip;
    if(cb_detector->HasElementFlags(ch, Detector_t::ElementFlag_t::BadTDC))
        return DoFitReturn_t::Skip;

    TH3* h_timewalk = dynamic_cast<TH3*>(hist);

    h_timewalk->GetZaxis()->SetRange(ch+1,ch+1);
    stringstream ss_name;
    ss_name << "Ch" << ch << "_yx";
    proj = dynamic_cast<TH2D*>(h_timewalk->Project3D("yx"));
    proj->FitSlicesY();
    means = dynamic_cast<TH1D*>(gDirectory->Get("timewalk_yx_1"));

    timewalks[ch]->Fit(means);

    last_timewalk = timewalks[ch]; // remember for display fit

    return DoFitReturn_t::Next;
}

void CB_TimeWalk::TheGUI::DisplayFit()
{
    c_fit->Show(means, last_timewalk.get(), true);

    c_extra->cd();
    proj->Draw("colz");
}

void CB_TimeWalk::TheGUI::StoreFit(unsigned channel)
{
    // the fit parameters contain the timewalk correction
    // and since we use pointers, the item in timewalks is already updated

    LOG(INFO) << "Stored Ch=" << channel;
}

bool CB_TimeWalk::TheGUI::FinishSlice()
{
    // don't request stop...
    return false;
}

void CB_TimeWalk::TheGUI::StoreFinishSlice(const interval<TID>& range)
{
    TCalibrationData cdata(
                GetName(),
                range.Start(),
                range.Stop()
                );

    // fill fit parameters
    for(unsigned ch=0;ch<cb_detector->GetNChannels();ch++) {
        const shared_ptr<gui::FitTimewalk>& func = timewalks[ch];
        cdata.FitParameters.emplace_back(ch, func->Save());
    }

    calibrationManager->Add(cdata, Calibration::AddMode_t::StrictRange);
}

