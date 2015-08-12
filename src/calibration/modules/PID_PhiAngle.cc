#include "PID_PhiAngle.h"

#include "DataManager.h"
#include "gui/CalCanvas.h"
#include "fitfunctions/FitGaus.h"

#include "expconfig/detectors/PID.h"
#include "base/Logger.h"
#include "base/std_ext.h"

#include "TGraph.h"

#include <limits>
#include <cmath>

using namespace ant;
using namespace ant::calibration;
using namespace ant::analysis;
using namespace ant::analysis::data;
using namespace std;

PID_PhiAngle::ThePhysics::ThePhysics(const string& name, unsigned nChannels) :
    Physics(name)
{
    const BinSettings pid_channels(nChannels);
    const BinSettings phibins(1000, -180, 3*180);

    pid_cb_phi_corr = HistFac.makeTH2D("CB/PID Cluster/Channel Correlation", "CB Cluster Phi / degree", "#",
                                       phibins, pid_channels, "pid_cb_phi_corr");
}

void PID_PhiAngle::ThePhysics::ProcessEvent(const Event& event)
{
    const auto& cands = event.Reconstructed().Candidates();

    // search for events with
    // one cluster in CB, one cluster in PID
    // ignore the matched candidates, since this is what
    // we want to calibrate

    const Cluster* cluster_pid = nullptr;
    double phi_cb = numeric_limits<double>::quiet_NaN();

    for(const auto& cand : cands) {

        auto cl_cb_  = cand->FindFirstCluster(Detector_t::Type_t::CB);

        if(cl_cb_ != nullptr) {
            // found more than one CB cluster
            if(isfinite(phi_cb))
                return;
            phi_cb = cand->Phi();
        }

        auto cl_pid_ = cand->FindFirstCluster(Detector_t::Type_t::PID);

        if(cl_pid_ != nullptr) {
            // found more than one pid cluster
            if(cluster_pid != nullptr)
                return;
            cluster_pid = cl_pid_;
        }
    }

    if(!isfinite(phi_cb) || cluster_pid == nullptr)
        return;

    phi_cb = std_ext::radian_to_degree(phi_cb);

    pid_cb_phi_corr->Fill(phi_cb,     cluster_pid->CentralElement);
    pid_cb_phi_corr->Fill(phi_cb+360, cluster_pid->CentralElement);

}

void PID_PhiAngle::ThePhysics::Finish()
{

}

void PID_PhiAngle::ThePhysics::ShowResult()
{
    canvas(GetName()) << drawoption("colz") << pid_cb_phi_corr << endc;
}

PID_PhiAngle::PID_PhiAngle(shared_ptr<expconfig::detector::PID> pid) :
    Module("PID_PhiAngle"),
    pid_detector(pid)
{
}

PID_PhiAngle::~PID_PhiAngle()
{
}

std::unique_ptr<analysis::Physics> PID_PhiAngle::GetPhysicsModule() {
    return std_ext::make_unique<ThePhysics>(GetName(), pid_detector->GetNChannels());
}

void PID_PhiAngle::GetGUIs(std::list<std::unique_ptr<gui::Manager_traits> >& guis) {
    guis.emplace_back(std_ext::make_unique<TheGUI>(GetName(), calibrationManager, pid_detector));
}


std::vector<std::list<TID> > ant::calibration::PID_PhiAngle::GetChangePoints() const
{
    return {calibrationManager->GetChangePoints(GetName())};
}

void ant::calibration::PID_PhiAngle::Update(size_t, const TID& id)
{
    TCalibrationData cdata;
    if(!calibrationManager->GetData(GetName(), id, cdata))
        return;
    if(cdata.Data.size() != 1)
        return;
    const TKeyValue<double>& kv = cdata.Data.front();
    pid_detector->SetPhiOffset(kv.Value);
}


PID_PhiAngle::TheGUI::TheGUI(const string& basename,
                             const std::shared_ptr<DataManager>& calmgr,
                             const std::shared_ptr<expconfig::detector::PID>& pid) :
    gui::Manager_traits(basename),
    calibrationManager(calmgr),
    pid_detector(pid)
{
}

string ant::calibration::PID_PhiAngle::TheGUI::GetHistogramName() const
{
    return GetName()+"/pid_cb_phi_corr";
}

unsigned ant::calibration::PID_PhiAngle::TheGUI::GetNumberOfChannels() const
{
    return pid_detector->GetNChannels();
}

void ant::calibration::PID_PhiAngle::TheGUI::InitGUI()
{
    c_singlechannel = new gui::CalCanvas("c_singlechannel", GetName()+": Single Channel");
    c_result = new gui::CalCanvas("c_result", GetName()+": Result");
}

std::list<gui::CalCanvas*> ant::calibration::PID_PhiAngle::TheGUI::GetCanvases() const
{
    return {c_singlechannel, c_result};
}

void ant::calibration::PID_PhiAngle::TheGUI::StartRange(const interval<TID>& range)
{
    // ask the detector for some reasonable starting values
    angles.resize(GetNumberOfChannels());
    for(size_t ch=0;ch<GetNumberOfChannels();ch++)
        angles[ch] = std_ext::radian_to_degree(pid_detector->GetPosition(ch).Phi());

    TCalibrationData cdata;
    if(calibrationManager->GetData(GetName()+"/SingleChannels", range.Start(), cdata)) {
        for(const TKeyValue<double>& kv : cdata.Data) {
            angles[kv.Key] = kv.Value;
        }
        for(const TKeyValue<vector<double>>& kv : cdata.FitParameters) {
            fitParameters.insert(make_pair(kv.Key, kv.Value));
        }
        LOG(INFO) << GetName() << ": Loaded previous single channel positions from database";
    }
    else {
        LOG(INFO) << GetName() << ": No previous data found";
    }

    // save a copy for comparison at finish stage
    previousAngles = angles;
}



gui::Manager_traits::DoFitReturn_t ant::calibration::PID_PhiAngle::TheGUI::DoFit(TH1* hist, unsigned channel)
{
    TH2* hist2 = dynamic_cast<TH2*>(hist);

    h_projection = hist2->ProjectionX("",channel,channel+1);

    func->SetDefaults(h_projection);
    const auto it_fit_param = fitParameters.find(channel);
    if(it_fit_param != fitParameters.end()) {
        VLOG(5) << "Loading previous fit parameters for channel " << channel;
        func->Load(it_fit_param->second);
    }

    func->Fit(h_projection);

    // always request display
    return DoFitReturn_t::Display;
}

void ant::calibration::PID_PhiAngle::TheGUI::DisplayFit()
{
    c_singlechannel->Show(h_projection, func.get());
}

void ant::calibration::PID_PhiAngle::TheGUI::StoreFit(unsigned channel)
{
    const double oldAngle = previousAngles[channel];

    double newAngle = func->GetPeakPosition();
    if(newAngle>180.0)
        newAngle -= 360.0;

    angles[channel] = newAngle;

    LOG(INFO) << "Stored Ch=" << channel << ": Angle " << newAngle << " from " << oldAngle;

    // don't forget the fit parameters
    fitParameters[channel] = func->Save();

    c_singlechannel->Clear();
    c_singlechannel->Update();
}

bool ant::calibration::PID_PhiAngle::TheGUI::FinishRange()
{
   h_result = new TGraph(GetNumberOfChannels());
   h_result->SetTitle("Result");
   h_result->GetXaxis()->SetTitle("Channel number");
   h_result->GetYaxis()->SetTitle("CB Phi position / degrees");
   for(size_t ch=0;ch<GetNumberOfChannels();ch++)
       h_result->SetPoint(ch, ch, angles[ch]);
   h_result->Fit("pol1","Q");

   c_result->cd();
   h_result->Draw();

   return true;
}

void ant::calibration::PID_PhiAngle::TheGUI::StoreFinishRange(const interval<TID>& range)
{
    delete h_result;
    c_result->Clear();
    c_result->Update();

    TCalibrationData cdata(
                "Unknown", /// \todo get static information about author/comment?
                "No Comment",
                time(nullptr),
                GetName()+"/SingleChannels",
                range.Start(),
                range.Stop()
                );

    // fill data
    for(unsigned ch=0;ch<angles.size();ch++) {
        cdata.Data.emplace_back(ch, angles[ch]);
    }

    // fill fit parameters (if any)
    for(const auto& it_map : fitParameters) {
        const unsigned ch = it_map.first;
        const vector<double>& params = it_map.second;
        cdata.FitParameters.emplace_back(ch, params);
    }

    calibrationManager->Add(cdata);

    LOG(INFO) << "Added TCalibrationData " << cdata;
}