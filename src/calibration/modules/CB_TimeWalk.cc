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
#include "TF1.h"
#include "base/std_ext/math.h"
#include "base/TH_ext.h"

#include <tree/TCluster.h>

#include <limits>
#include <cmath>

using namespace ant;
using namespace ant::calibration;
using namespace std;

CB_TimeWalk::CB_TimeWalk(
        const shared_ptr<expconfig::detector::CB>& cb,
        const shared_ptr<DataManager>& calmgr,
        const interval<double>& timeWindow,
        double badTDC_EnergyThreshold
        ) :
    Module("CB_TimeWalk"),
    cb_detector(cb),
    calibrationManager(calmgr),
    TimeWindow(timeWindow),
    BadTDC_EnergyThreshold(badTDC_EnergyThreshold)
{
    for(unsigned ch=0;ch<cb_detector->GetNChannels();ch++) {
        timewalks.emplace_back(make_shared<gui::FitTimewalk>());
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


    for(TClusterHit& clusterhit : clusterhits) {

        // check if this hit belongs to a bad TDC chnanel
        if(cb_detector->HasElementFlags(clusterhit.Channel, Detector_t::ElementFlag_t::BadTDC)) {
            // if energy is above threhshold, keep the clusterhit
            // otherwise invalidate its energy (then ignored by clustering)
            if(clusterhit.Energy < BadTDC_EnergyThreshold)
                clusterhit.Energy = std_ext::NaN;
            continue;
        }

        // do timewalk correction
        // use uncalibrated energy for that
        // to stay independent of energy calibration


        double raw_energy = std_ext::NaN;
        vector<double> timings;
        for(auto& datum : clusterhit.Data) {
            if(datum.Type == Channel_t::Type_t::Integral) {
                raw_energy = datum.Value.Uncalibrated;
            }
            else if(datum.Type == Channel_t::Type_t::Timing) {
                timings.push_back(datum.Value.Calibrated);
            }
        }

        // don't bother with stuff where no enerygy or time info is present
        if(!isfinite(raw_energy)) {
            VLOG(7) << "Found " << clusterhit << " without raw energy information.";
            continue;
        }

        if(timings.empty()) {
            VLOG(7) << "Found " << clusterhit << " without any timings.";
            continue;
        }

        // Eval of Timewalk function handles
        // case of low raw energy for us
        auto deltaT = timewalks[clusterhit.Channel]->Eval(std::log10(raw_energy));

        // still check if we got a finite value
        // else just do nothing with clusterhit
        if(isfinite(deltaT)) {

            // find the timing which is closest to deltaT
            auto it_min_timing = std::min_element(
                                     timings.begin(), timings.end(),
                                     [deltaT] (double a, double b) {
                return abs(a - deltaT) < abs(b - deltaT);
            });

            clusterhit.Time = *it_min_timing - deltaT;

            // get rid of clusterhit if outside timewindow
            if(!TimeWindow.Contains(clusterhit.Time)) {
                clusterhit.Time = std_ext::NaN;
            }
        }
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
    slicesY_gaus = new TF1("slicesY_gaus","gaus");
}

shared_ptr<TH1> CB_TimeWalk::TheGUI::GetHistogram(const WrapTFile& file) const
{
    return file.GetSharedHist<TH1>(GetName()+"/timewalk");
}

unsigned CB_TimeWalk::TheGUI::GetNumberOfChannels() const
{
    return cb_detector->GetNChannels();
}

void CB_TimeWalk::TheGUI::InitGUI(gui::ManagerWindow_traits& window)
{
    c_fit = window.AddCalCanvas();
    c_extra = window.AddCalCanvas();

    window.AddNumberEntry("Chi2/NDF limit for autostop", AutoStopOnChi2);
    window.AddNumberEntry("SlicesYEntryCut", slicesY_entryCut);
    window.AddNumberEntry("SlicesYIQRFactor low  (outlier detection)", slicesY_IQRFactor_lo);
    window.AddNumberEntry("SlicesYIQRFactor high (outlier detection)", slicesY_IQRFactor_hi);
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
            LOG(ERROR) << "Ignoring too large channel key=" << kv.Key;
            continue;
        }
        fitParameters[kv.Key] = kv.Value;
    }
}


gui::CalibModule_traits::DoFitReturn_t CB_TimeWalk::TheGUI::DoFit(const TH1& hist, unsigned ch)
{
    if(cb_detector->IsIgnored(ch))
        return DoFitReturn_t::Skip;
    if(cb_detector->HasElementFlags(ch, Detector_t::ElementFlag_t::BadTDC))
        return DoFitReturn_t::Skip;

    auto& h_timewalk = dynamic_cast<const TH3&>(hist);

    // cannot call SetRange on const TAxis
    TH3* htmp = (TH3*) h_timewalk.Clone(TString::Format("%s_clone", h_timewalk.GetName()));
    htmp->GetZaxis()->SetRange(ch+1,ch+1);
    proj = dynamic_cast<TH2D*>(htmp->Project3D("yx"));
    delete htmp;

    means = TH_ext::FitSlicesY(proj, slicesY_gaus, slicesY_entryCut,
                               slicesY_IQRFactor_lo, slicesY_IQRFactor_hi);
    means->SetMinimum(proj->GetYaxis()->GetXmin());
    means->SetMaximum(proj->GetYaxis()->GetXmax());

    auto& func = timewalks[ch];
    func->SetDefaults(means);
    func->SetRange({1, means->GetXaxis()->GetXmax()});
    const auto it_fit_param = fitParameters.find(ch);
    if(it_fit_param != fitParameters.end()) {
        VLOG(5) << "Loading previous fit parameters for channel " << ch;
        func->Load(it_fit_param->second);
    }
    last_timewalk = func; // remember for display fit


    auto fit_loop = [this,func] (size_t retries) {
        do {
            func->Fit(means);
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

void CB_TimeWalk::TheGUI::DisplayFit()
{
    c_fit->Show(means, last_timewalk.get(), true);

    c_extra->cd();
    proj->Draw("colz");
    /// \todo make func update when other canvas is edited,
    /// clicking on it is enough to redraw function
    last_timewalk->Draw();
}

void CB_TimeWalk::TheGUI::StoreFit(unsigned channel)
{
    auto& func = timewalks[channel];
    fitParameters[channel] = func->Save();
    // the fit parameters contain the timewalk correction
    // and since we use pointers, the item in timewalks is already updated

    LOG(INFO) << "Stored Ch=" << channel << " Parameters: " << fitParameters[channel];
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

    calibrationManager->Add(cdata);
}

