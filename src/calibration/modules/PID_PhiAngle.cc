#include "PID_PhiAngle.h"

#include "DataManager.h"
#include "gui/CalCanvas.h"
#include "fitfunctions/FitGaus.h"
#include "tree/TCalibrationData.h"

#include "expconfig/detectors/PID.h"
#include "base/Logger.h"
#include "base/std_ext/math.h"

#include "TGraph.h"
#include "TFitResult.h"

#include <limits>
#include <cmath>

#include "TF1.h"

using namespace ant;
using namespace ant::calibration;
using namespace ant::analysis;
using namespace ant::analysis::data;
using namespace std;

PID_PhiAngle::ThePhysics::ThePhysics(const string& name, unsigned nChannels) :
    Physics(name),
    theta_range(40.0*TMath::DegToRad(), 140*TMath::DegToRad())
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

        if(!theta_range.Contains(cand->Theta()))
            continue;

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

    /// \todo search all clusters, leave candidates alone
    for(const Cluster& cl : event.Reconstructed().AllClusters()) {
        if(cl.Detector != Detector_t::Type_t::PID)
            continue;
        if(!isfinite(cl.Energy) || !isfinite(cl.Time))
            continue;
        // found more than one PID cluster
        if(cluster_pid != nullptr)
            return;
        cluster_pid = addressof(cl);
    }

    if(!isfinite(phi_cb) || cluster_pid == nullptr)
        return;

    const double phi_cb_degrees = std_ext::radian_to_degree(phi_cb);

    pid_cb_phi_corr->Fill(phi_cb_degrees,     cluster_pid->CentralElement);
    pid_cb_phi_corr->Fill(phi_cb_degrees+360, cluster_pid->CentralElement);

}

void PID_PhiAngle::ThePhysics::Finish()
{

}

void PID_PhiAngle::ThePhysics::ShowResult()
{
    canvas(GetName()) << drawoption("colz") << pid_cb_phi_corr << endc;
}

PID_PhiAngle::PID_PhiAngle(
        const std::shared_ptr<expconfig::detector::PID>& pid,
        const std::shared_ptr<DataManager>& calmgr) :
    Module("PID_PhiAngle"),
    pid_detector(pid),
    calibrationManager(calmgr)
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


std::vector<std::list<TID> > PID_PhiAngle::GetChangePoints() const
{
    return {calibrationManager->GetChangePoints(GetName())};
}

void PID_PhiAngle::Update(size_t, const TID& id)
{
    TCalibrationData cdata;
    if(!calibrationManager->GetData(GetName(), id, cdata))
        return;
    if(cdata.Data.size() != 1)
        return;
    const TKeyValue<double>& kv = cdata.Data.front();
    pid_detector->SetPhiOffset(kv.Value);
}


/**
 * @brief The PID_PhiAngle::TheGUI::_FitGauss: override the SetDefaults for PID phi angle fits
 */
class PID_PhiAngle::TheGUI::_FitGauss: public gui::FitGaus {
public:
    using gui::FitGaus::FitGaus;

    static constexpr const auto halfrange = 50.0;

    virtual void SetDefaults(TH1 *hist) override {
        if(hist) {

            func->SetParameter(0,hist->GetMaximum());

            double max_pos = hist->GetXaxis()->GetBinCenter(hist->GetMaximumBin());

            if(max_pos < -180+halfrange)
                max_pos += 360.0;

            if(max_pos > 540-halfrange)
                max_pos -= 360.0;

            func->SetParameter(1,max_pos);
            SetRange({max_pos-halfrange,max_pos+halfrange});

         } else {
            SetRange({0,200});
            func->SetParameter(0,0.8);
            func->SetParameter(1,100);
        }
        func->SetParameter(2,10.0);
    }
};


PID_PhiAngle::TheGUI::TheGUI(const string& basename,
                             const std::shared_ptr<DataManager>& calmgr,
                             const std::shared_ptr<expconfig::detector::PID>& pid) :
    gui::Manager_traits(basename),
    calibrationManager(calmgr),
    pid_detector(pid),
    func(make_shared<_FitGauss>())
{
}

PID_PhiAngle::TheGUI::~TheGUI()
{

}

string PID_PhiAngle::TheGUI::GetHistogramName() const
{
    return GetName()+"/pid_cb_phi_corr";
}

unsigned PID_PhiAngle::TheGUI::GetNumberOfChannels() const
{
    return pid_detector->GetNChannels();
}

void PID_PhiAngle::TheGUI::InitGUI(gui::ManagerWindow_traits* window)
{
    canvas = window->AddCalCanvas();
}

void PID_PhiAngle::TheGUI::StartRange(const interval<TID>& range)
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



gui::Manager_traits::DoFitReturn_t PID_PhiAngle::TheGUI::DoFit(TH1* hist, unsigned channel,
                                                               const DoFitOptions_t& options)
{
    TH2* hist2 = dynamic_cast<TH2*>(hist);

    h_projection = hist2->ProjectionX("",channel+1,channel+1);

    func->SetDefaults(h_projection);
    const auto it_fit_param = fitParameters.find(channel);
    if(it_fit_param != fitParameters.end()
       && !options.IgnorePreviousFitParameters) {
        VLOG(5) << "Loading previous fit parameters for channel " << channel;
        func->Load(it_fit_param->second);
    }

    func->Fit(h_projection);

    // always request display
    return DoFitReturn_t::Next;
}

void PID_PhiAngle::TheGUI::DisplayFit()
{
    canvas->Show(h_projection, func.get());
}

void PID_PhiAngle::TheGUI::StoreFit(unsigned channel)
{
    const double oldAngle = previousAngles[channel];

    double newAngle = func->GetPeakPosition();
    if(newAngle>180.0)
        newAngle -= 360.0;
    if(newAngle<0)
        newAngle += 360.0;


    angles[channel] = newAngle;

    LOG(INFO) << "Stored Ch=" << channel << ": Angle " << newAngle << " from " << oldAngle;

    // don't forget the fit parameters
    fitParameters[channel] = func->Save();
}

//double Pol1Wrap360(const double* x, const double* p) {
//    const auto& offset = p[0];
//    const auto& slope  = p[1];

//    auto y = offset + slope*x[0];

//    const int n = floor(y / 360.0);
//    y -= n * 360.0;

//    return y;
//}

bool PID_PhiAngle::TheGUI::FinishRange()
{
   h_result = new TGraph(GetNumberOfChannels());

   for(size_t ch=0;ch<GetNumberOfChannels();ch++) {

       auto angle = angles[ch];

       // always make sure values increase (wrap around)
       /// @todo what to do if slope is negative?
       if(ch!=0)
           while(angles[ch-1] > angle)
               angle += 360.0;

       h_result->SetPoint(ch, ch, angle);
   }

   TFitResultPtr r = h_result->Fit("pol1","QS");
   phi_offset = r->Value(0);
   LOG(INFO) << "Found Phi offset of first channel: " << phi_offset;

   h_result->SetTitle("Result");
   h_result->GetXaxis()->SetTitle("Channel number");
   h_result->GetYaxis()->SetTitle("CB Phi position / degrees");
   h_result->SetMarkerSize(2);
   h_result->SetMarkerStyle(kMultiply);

   canvas->cd();
   h_result->Draw("AP");

   return true;
}

void PID_PhiAngle::TheGUI::StoreFinishRange(const interval<TID>& range)
{
    delete h_result;

    TCalibrationData cdata(
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

    TCalibrationData cdata_offset(
                GetName(),
                range.Start(),
                range.Stop()
                );
    cdata_offset.Data.emplace_back(0, phi_offset);


    calibrationManager->Add(cdata);
    calibrationManager->Add(cdata_offset);
}

