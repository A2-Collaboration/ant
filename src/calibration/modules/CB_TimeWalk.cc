#include "CB_TimeWalk.h"

#include "DataManager.h"
#include "gui/CalCanvas.h"
#include "fitfunctions/FitTimewalk.h"
#include "tree/TCalibrationData.h"

#include "reconstruct/Clustering.h"

#include "expconfig/detectors/CB.h"
#include "base/Logger.h"

#include "TGraph.h"
#include "TFitResult.h"
#include "TObjArray.h"
#include "TDirectory.h"

#include <tree/TCluster.h>

#include <limits>
#include <cmath>

using namespace ant;
using namespace ant::calibration;
using namespace ant::analysis;
using namespace ant::analysis::data;
using namespace std;

CB_TimeWalk::ThePhysics::ThePhysics(const string& name, const std::shared_ptr<expconfig::detector::CB>& cb) :
    Physics(name),
    cb_detector(cb)
{
    BinSettings bins_energy(400,0,500);
    BinSettings bins_channels(cb_detector->GetNChannels());

    h_timewalk =
            HistFac.makeTH3D(
                "CB TimeWalk",
                "Energy / MeV",
                "Time / ns",
                "Channel",
                bins_energy,
                BinSettings(100,-100,100),
                bins_channels,
                "timewalk"
                );
    h_timewalk_overview =
            HistFac.makeTH2D(
                "CB Timewalk Overview",
                "Energy / MeV",
                "Channel",
                bins_energy,
                bins_channels,
                "timewalk_overview"
                );
}

void CB_TimeWalk::ThePhysics::ProcessEvent(const Event& event)
{
    for(const auto& cand: event.Reconstructed().Candidates()) {
        for(const Cluster& cluster: cand->Clusters) {
            if(cluster.Detector != Detector_t::Type_t::CB)
                continue;
            for(const Cluster::Hit& hit : cluster.Hits) {
                // found the hit of the central element
                // now search for its timing information
                double time = numeric_limits<double>::quiet_NaN();
                double energy = numeric_limits<double>::quiet_NaN();
                for(const Cluster::Hit::Datum& d : hit.Data) {
                    if(d.Type == Channel_t::Type_t::Timing)
                        time = d.Value;
                    if(d.Type == Channel_t::Type_t::Integral)
                        energy = d.Value;
                }
                h_timewalk->Fill(energy, time, hit.Channel);
            }
        }
    }
}

void CB_TimeWalk::ThePhysics::Finish()
{

}

void CB_TimeWalk::ThePhysics::ShowResult()
{
    for(unsigned ch=0;ch<cb_detector->GetNChannels();ch++) {
        if(cb_detector->IsIgnored(ch))
            continue;
        LOG(INFO) << "Fitting Channel=" << ch;
        h_timewalk->GetZaxis()->SetRange(ch,ch+1);
        stringstream ss_name;
        ss_name << "Ch" << ch << "_yx";
        TH2* proj = dynamic_cast<TH2*>(h_timewalk->Project3D(ss_name.str().c_str()));
        TObjArray aSlices;
        proj->FitSlicesY(nullptr, 0, -1, 0, "QNR", &aSlices);
        TH1D* means = dynamic_cast<TH1D*>(aSlices.At(1));
        for(Int_t x=0;x<means->GetNbinsX()+1;x++) {
            h_timewalk_overview->SetBinContent(x, ch+1, means->GetBinContent(x));
        }
        delete proj;
    }
    canvas(GetName()) << drawoption("colz") << h_timewalk_overview << endc;
}

CB_TimeWalk::CB_TimeWalk(const shared_ptr<expconfig::detector::CB>& cb,
        const shared_ptr<DataManager>& calmgr,
                         const interval<double>& timeWindow) :
    Module("CB_TimeWalk"),
    cb_detector(cb),
    calibrationManager(calmgr),
    TimeWindow(timeWindow)
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

    list< reconstruct::AdaptorTClusterHit >& clusterhits = it_sorted_clusterhits->second;

    auto it_clusterhit = clusterhits.begin();

    while(it_clusterhit != clusterhits.end()) {
        reconstruct::AdaptorTClusterHit& clusterhit = *it_clusterhit;
        // do timewalk correction
        clusterhit.Time -= timewalks[clusterhit.Hit->Channel]->Eval(clusterhit.Energy);
        // get rid of clusterhit if outside timewindow
        if(std::isfinite(clusterhit.Time) && !TimeWindow.Contains(clusterhit.Time))
            it_clusterhit = clusterhits.erase(it_clusterhit);
        else
            ++it_clusterhit;
    }
}

unique_ptr<analysis::Physics> CB_TimeWalk::GetPhysicsModule() {
    return std_ext::make_unique<ThePhysics>(GetName(), cb_detector);
}

void CB_TimeWalk::GetGUIs(list<unique_ptr<gui::Manager_traits> >& guis) {
    guis.emplace_back(std_ext::make_unique<TheGUI>(GetName(), calibrationManager, cb_detector, timewalks));
}


vector<list<TID> > CB_TimeWalk::GetChangePoints() const
{
    return {calibrationManager->GetChangePoints(GetName())};
}

void CB_TimeWalk::Update(size_t, const TID& id)
{
    TCalibrationData cdata;
    if(!calibrationManager->GetData(GetName(), id, cdata))
        return;
    for(const TKeyValue<vector<double>>& kv : cdata.FitParameters) {
        if(kv.Key>=timewalks.size()) {
            LOG(ERROR) << "Ignoring too large key=" << kv.Key;
            continue;
        }
        timewalks[kv.Key]->Load(kv.Value);
    }
}

void CB_TimeWalk::UpdatedTIDFlags(const TID& id)
{
    IsMC = id.isSet(TID::Flags_t::MC);
}


CB_TimeWalk::TheGUI::TheGUI(const string& basename,
                            const shared_ptr<DataManager>& calmgr,
                            const shared_ptr<expconfig::detector::CB>& cb,
                            std::vector< std::shared_ptr<gui::FitTimewalk> >& timewalks_) :
    gui::Manager_traits(basename),
    calibrationManager(calmgr),
    cb_detector(cb),
    timewalks(timewalks_)
{
}

string CB_TimeWalk::TheGUI::GetHistogramName() const
{
    return GetName()+"/timewalk";
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

void CB_TimeWalk::TheGUI::StartRange(const interval<TID>& range)
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



gui::Manager_traits::DoFitReturn_t CB_TimeWalk::TheGUI::DoFit(TH1* hist, unsigned ch,
                                                              const DoFitOptions_t&)
{
    if(cb_detector->IsIgnored(ch))
        return DoFitReturn_t::Skip;

    TH3* h_timewalk = dynamic_cast<TH3*>(hist);

    h_timewalk->GetZaxis()->SetRange(ch,ch+1);
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

bool CB_TimeWalk::TheGUI::FinishRange()
{
    // don't request stop...
    return false;
}

void CB_TimeWalk::TheGUI::StoreFinishRange(const interval<TID>& range)
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

