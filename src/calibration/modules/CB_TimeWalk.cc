#include "CB_TimeWalk.h"

#include "DataManager.h"
#include "gui/CalCanvas.h"
#include "fitfunctions/FitGaus.h"

#include "expconfig/detectors/CB.h"
#include "base/Logger.h"
#include "base/std_ext.h"

#include "TGraph.h"
#include "TFitResult.h"

#include <limits>
#include <cmath>

using namespace ant;
using namespace ant::calibration;
using namespace ant::analysis;
using namespace ant::analysis::data;
using namespace std;

CB_TimeWalk::ThePhysics::ThePhysics(const string& name, unsigned nChannels) :
    Physics(name)
{
    h_timewalk = HistFac.makeTH3D(
                     "CB TimeWalk",
                     "Energy / MeV",
                     "Time / ns",
                     "Channel",
                     BinSettings(1000,0,1000),
                     BinSettings(500,-100,100),
                     BinSettings(nChannels),
                     "timewalk"
                     );
}

void CB_TimeWalk::ThePhysics::ProcessEvent(const Event& event)
{
    for(const auto& cand: event.Reconstructed().Candidates()) {
        for(const Cluster& cluster: cand->Clusters) {
            if(cluster.Detector != Detector_t::Type_t::CB)
                continue;
            for(const Cluster::Hit& hit : cluster.Hits) {
                if(cluster.CentralElement != hit.Channel)
                    continue;
                // found the hit of the central element
                double time = numeric_limits<double>::quiet_NaN();
                double energy = numeric_limits<double>::quiet_NaN();
                for(const Cluster::Hit::Datum& d : hit.Data) {
                    if(d.Type == Channel_t::Type_t::Timing)
                        time = d.Value;
                    if(d.Type == Channel_t::Type_t::Integral)
                        energy = d.Value;
                }
                h_timewalk->Fill(energy, time, hit.Channel);
                return;
            }
        }
    }
}

void CB_TimeWalk::ThePhysics::Finish()
{

}

void CB_TimeWalk::ThePhysics::ShowResult()
{
    canvas(GetName()) << drawoption("colz") << h_timewalk << endc;
}

CB_TimeWalk::CB_TimeWalk(
        const std::shared_ptr<expconfig::detector::CB>& cb,
        const std::shared_ptr<DataManager>& calmgr) :
    Module("CB_TimeWalk"),
    cb_detector(cb),
    calibrationManager(calmgr)
{
}

CB_TimeWalk::~CB_TimeWalk()
{
}

std::unique_ptr<analysis::Physics> CB_TimeWalk::GetPhysicsModule() {
    return std_ext::make_unique<ThePhysics>(GetName(), cb_detector->GetNChannels());
}

void CB_TimeWalk::GetGUIs(std::list<std::unique_ptr<gui::Manager_traits> >& guis) {
    guis.emplace_back(std_ext::make_unique<TheGUI>(GetName(), calibrationManager, cb_detector));
}


std::vector<std::list<TID> > ant::calibration::CB_TimeWalk::GetChangePoints() const
{
    return {calibrationManager->GetChangePoints(GetName())};
}

void ant::calibration::CB_TimeWalk::Update(size_t, const TID& id)
{
//    TCalibrationData cdata;
//    if(!calibrationManager->GetData(GetName(), id, cdata))
//        return;
//    if(cdata.Data.size() != 1)
//        return;
//    const TKeyValue<double>& kv = cdata.Data.front();
//    pid_detector->SetPhiOffset(kv.Value);
}


CB_TimeWalk::TheGUI::TheGUI(const string& basename,
                             const std::shared_ptr<DataManager>& calmgr,
                             const std::shared_ptr<expconfig::detector::CB>& cb) :
    gui::Manager_traits(basename),
    calibrationManager(calmgr),
    cb_detector(cb),
    func(make_shared<gui::FitGaus>())
{
}

string ant::calibration::CB_TimeWalk::TheGUI::GetHistogramName() const
{
    return GetName()+"/pid_cb_phi_corr";
}

unsigned ant::calibration::CB_TimeWalk::TheGUI::GetNumberOfChannels() const
{
    return cb_detector->GetNChannels();
}

void ant::calibration::CB_TimeWalk::TheGUI::InitGUI()
{
//    c_singlechannel = new gui::CalCanvas("c_singlechannel", GetName()+": Single Channel");
//    c_result = new gui::CalCanvas("c_result", GetName()+": Result");
}

std::list<gui::CalCanvas*> ant::calibration::CB_TimeWalk::TheGUI::GetCanvases() const
{
    //return {c_singlechannel, c_result};
    return {};
}

void ant::calibration::CB_TimeWalk::TheGUI::StartRange(const interval<TID>& range)
{
//    // ask the detector for some reasonable starting values
//    angles.resize(GetNumberOfChannels());
//    for(size_t ch=0;ch<GetNumberOfChannels();ch++)
//        angles[ch] = std_ext::radian_to_degree(pid_detector->GetPosition(ch).Phi());

//    TCalibrationData cdata;
//    if(calibrationManager->GetData(GetName()+"/SingleChannels", range.Start(), cdata)) {
//        for(const TKeyValue<double>& kv : cdata.Data) {
//            angles[kv.Key] = kv.Value;
//        }
//        for(const TKeyValue<vector<double>>& kv : cdata.FitParameters) {
//            fitParameters.insert(make_pair(kv.Key, kv.Value));
//        }
//        LOG(INFO) << GetName() << ": Loaded previous single channel positions from database";
//    }
//    else {
//        LOG(INFO) << GetName() << ": No previous data found";
//    }

//    // save a copy for comparison at finish stage
//    previousAngles = angles;
}



gui::Manager_traits::DoFitReturn_t ant::calibration::CB_TimeWalk::TheGUI::DoFit(TH1* hist, unsigned channel)
{
//    TH2* hist2 = dynamic_cast<TH2*>(hist);

//    h_projection = hist2->ProjectionX("",channel,channel+1);

//    func->SetDefaults(h_projection);
//    const auto it_fit_param = fitParameters.find(channel);
//    if(it_fit_param != fitParameters.end()) {
//        VLOG(5) << "Loading previous fit parameters for channel " << channel;
//        func->Load(it_fit_param->second);
//    }

//    func->Fit(h_projection);

    // always request display
    return DoFitReturn_t::Display;
}

void ant::calibration::CB_TimeWalk::TheGUI::DisplayFit()
{
    //c_singlechannel->Show(h_projection, func.get());
}

void ant::calibration::CB_TimeWalk::TheGUI::StoreFit(unsigned channel)
{
//    const double oldAngle = previousAngles[channel];

//    double newAngle = func->GetPeakPosition();
//    if(newAngle>180.0)
//        newAngle -= 360.0;
//    if(newAngle<0)
//        newAngle += 360.0;


//    angles[channel] = newAngle;

//    LOG(INFO) << "Stored Ch=" << channel << ": Angle " << newAngle << " from " << oldAngle;

//    // don't forget the fit parameters
//    fitParameters[channel] = func->Save();

//    c_singlechannel->Clear();
//    c_singlechannel->Update();
}

bool ant::calibration::CB_TimeWalk::TheGUI::FinishRange()
{
//   h_result = new TGraph(GetNumberOfChannels());
//   for(size_t ch=0;ch<GetNumberOfChannels();ch++)
//       h_result->SetPoint(ch, ch, angles[ch]);
//   TFitResultPtr r = h_result->Fit("pol1","QS");
//   phi_offset = r->Value(0);
//   LOG(INFO) << "Found Phi offset of first channel: " << phi_offset;

//   h_result->SetTitle("Result");
//   h_result->GetXaxis()->SetTitle("Channel number");
//   h_result->GetYaxis()->SetTitle("CB Phi position / degrees");
//   h_result->SetMarkerSize(2);
//   h_result->SetMarkerStyle(kMultiply);

//   c_result->cd();
//   h_result->Draw("AP");

//   c_result->Modified();
//   c_result->Update();

   return true;
}

void ant::calibration::CB_TimeWalk::TheGUI::StoreFinishRange(const interval<TID>& range)
{
//    delete h_result;
//    c_result->Clear();
//    c_result->Update();

//    TCalibrationData cdata(
//                "Unknown", /// \todo get static information about author/comment?
//                "No Comment",
//                time(nullptr),
//                GetName()+"/SingleChannels",
//                range.Start(),
//                range.Stop()
//                );

//    // fill data
//    for(unsigned ch=0;ch<angles.size();ch++) {
//        cdata.Data.emplace_back(ch, angles[ch]);
//    }

//    // fill fit parameters (if any)
//    for(const auto& it_map : fitParameters) {
//        const unsigned ch = it_map.first;
//        const vector<double>& params = it_map.second;
//        cdata.FitParameters.emplace_back(ch, params);
//    }

//    TCalibrationData cdata_offset(
//                "Unknown", /// \todo get static information about author/comment?
//                "No Comment",
//                time(nullptr),
//                GetName(),
//                range.Start(),
//                range.Stop()
//                );
//    cdata_offset.Data.emplace_back(0, phi_offset);


//    calibrationManager->Add(cdata);
//    calibrationManager->Add(cdata_offset);
//    LOG(INFO) << "Added TCalibrationData: " << cdata;
//    LOG(INFO) << "Added TCalibrationData: " << cdata_offset;
}