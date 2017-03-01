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

void CB_TimeWalk::TheGUI::InitGUI(gui::ManagerWindow_traits* window)
{
    c_fit = window->AddCalCanvas();
    c_extra = window->AddCalCanvas();
    c_fit->ConnectOtherCalCanvas(c_extra); // update c_extra if c_fit is updated

    window->AddNumberEntry("Chi2/NDF limit for autostop", AutoStopOnChi2);
    window->AddNumberEntry("SlicesYEntryCut", slicesY_entryCut);
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
        timewalks[kv.Key]->Load(kv.Value);
    }
}

// copied and adapted from TH2::FitSlicesY/DoFitSlices
void MyFitSlicesY(TH2* h, TF1 *f1, Int_t cut)
{
    TAxis& outerAxis = *h->GetXaxis();
    Int_t nbins  = outerAxis.GetNbins();
    Int_t firstbin = 0;
    Int_t lastbin = nbins + 1;
    Int_t ngroup = 1;
    Int_t nstep = ngroup;

    Int_t npar = f1->GetNpar();
    if (npar <= 0) return;

    //Create one histogram for each function parameter
    Int_t ipar;
    TH1D **hlist = new TH1D*[npar];
    char *name   = new char[2000];
    char *title  = new char[2000];
    const TArrayD *bins = outerAxis.GetXbins();
    for (ipar=0;ipar<npar;ipar++) {
        snprintf(name,2000,"%s_%d",h->GetName(),ipar);
        snprintf(title,2000,"Fitted value of par[%d]=%s",ipar,f1->GetParName(ipar));
        delete gDirectory->FindObject(name);
        if (bins->fN == 0) {
            hlist[ipar] = new TH1D(name,title, nbins, outerAxis.GetXmin(), outerAxis.GetXmax());
        } else {
            hlist[ipar] = new TH1D(name,title, nbins,bins->fArray);
        }
        hlist[ipar]->GetXaxis()->SetTitle(outerAxis.GetTitle());
    }

    //Loop on all bins in Y, generate a projection along X
    Int_t bin;
    Long64_t nentries;
    // in case of sliding merge nstep=1, i.e. do slices starting for every bin
    // now do not slices case with overflow (makes more sense)
    for (bin=firstbin;bin+ngroup-1<=lastbin;bin += nstep) {
        TH1D *hp = h->ProjectionY("_temp",bin,bin+ngroup-1,"e");
        if (hp == 0) continue;
        nentries = Long64_t(hp->GetEntries());
        if (nentries == 0 || nentries < cut) {delete hp; continue;}


        const double max_pos = hp->GetXaxis()->GetBinCenter(hp->GetMaximumBin());
        const double sigma = hp->GetRMS();

        // setting meaningful start parameters and limits
        // is crucial for a good fit!
        f1->SetParameter(0, hp->GetMaximum());
        f1->SetParLimits(0, 0, hp->GetMaximum());
        f1->SetParLimits(1, max_pos-4*sigma, max_pos+4*sigma);
        f1->SetParameter(1, max_pos);
        f1->SetParLimits(2, 0, 60);
        f1->SetParameter(2, sigma); // set sigma
        f1->SetRange(max_pos-4*sigma, max_pos+4*sigma);

        hp->Fit(f1,"QBNR"); // B important for limits!



        Int_t npfits = f1->GetNumberFitPoints();
        if (npfits > npar && npfits >= cut) {
            if(f1->GetParError(2)<10) {

                Int_t binOn = bin + ngroup/2;
                for (ipar=0;ipar<npar;ipar++) {
                    hlist[ipar]->Fill(outerAxis.GetBinCenter(binOn),f1->GetParameter(ipar));
                    hlist[ipar]->SetBinError(binOn,f1->GetParError(ipar));
                }
            }
        }
        delete hp;
    }
    delete [] name;
    delete [] title;
    delete [] hlist;
}


gui::CalibModule_traits::DoFitReturn_t CB_TimeWalk::TheGUI::DoFit(TH1* hist, unsigned ch)
{
    if(cb_detector->IsIgnored(ch))
        return DoFitReturn_t::Skip;
    if(cb_detector->HasElementFlags(ch, Detector_t::ElementFlag_t::BadTDC))
        return DoFitReturn_t::Skip;

    TH3* h_timewalk = dynamic_cast<TH3*>(hist);

    h_timewalk->GetZaxis()->SetRange(ch+1,ch+1);
    proj = dynamic_cast<TH2D*>(h_timewalk->Project3D("yx"));

    MyFitSlicesY(proj, slicesY_gaus, slicesY_entryCut);
    means = dynamic_cast<TH1D*>(gDirectory->Get("timewalk_yx_1"));

    means->SetMinimum(proj->GetYaxis()->GetXmin());
    means->SetMaximum(proj->GetYaxis()->GetXmax());
    auto& func = timewalks[ch];
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
    last_timewalk->Draw();
}

void CB_TimeWalk::TheGUI::StoreFit(unsigned channel)
{
    // the fit parameters contain the timewalk correction
    // and since we use pointers, the item in timewalks is already updated

    LOG(INFO) << "Stored Ch=" << channel << " Parameters: " << timewalks[channel]->Save();
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

