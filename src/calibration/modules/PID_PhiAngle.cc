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
#include "TH1.h"
#include "TH2.h"

using namespace std;
using namespace ant;
using namespace ant::calibration;

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

void PID_PhiAngle::GetGUIs(std::list<std::unique_ptr<gui::CalibModule_traits> >& guis, OptionsPtr) {
    guis.emplace_back(std_ext::make_unique<TheGUI>(GetName(), calibrationManager, pid_detector));
}

std::list<Updateable_traits::Loader_t> PID_PhiAngle::GetLoaders()
{
    return {
        [this] (const TID& currPoint, TID& nextChangePoint) {
            TCalibrationData cdata;
            if(!calibrationManager->GetData(GetName(), currPoint, cdata, nextChangePoint))
                return;
            if(cdata.Data.size() != 1)
                return;
            const TKeyValue<double>& kv = cdata.Data.front();
            pid_detector->SetPhiOffset(kv.Value);
        }
    };
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
    gui::CalibModule_traits(basename),
    calibrationManager(calmgr),
    pid_detector(pid),
    func(make_shared<_FitGauss>())
{
}

PID_PhiAngle::TheGUI::~TheGUI()
{

}

shared_ptr<TH1> PID_PhiAngle::TheGUI::GetHistogram(const WrapTFile& file) const
{
    return file.GetSharedHist<TH1>(GetName()+"/pid_cb_phi_corr");
}

unsigned PID_PhiAngle::TheGUI::GetNumberOfChannels() const
{
    return pid_detector->GetNChannels();
}

void PID_PhiAngle::TheGUI::InitGUI(gui::ManagerWindow_traits& window)
{
    window.AddCheckBox("Ignore prev fit params", IgnorePreviousFitParameters);
    canvas = window.AddCalCanvas();
}

void PID_PhiAngle::TheGUI::StartSlice(const interval<TID>& range)
{
    // ask the detector for some reasonable starting values
    angles.resize(GetNumberOfChannels());
    for(size_t ch=0;ch<GetNumberOfChannels();ch++)
        angles[ch] = std_ext::radian_to_degree(pid_detector->GetPosition(ch).Phi());

    TCalibrationData cdata;
    if(calibrationManager->GetData(GetName()+"_SingleChannels", range.Start(), cdata)) {
        for(const TKeyValue<double>& kv : cdata.Data) {
            angles[kv.Key] = kv.Value;
        }
        for(const TKeyValue<vector<double>>& kv : cdata.FitParameters) {
            fitParameters[kv.Key] = kv.Value;
        }
        LOG(INFO) << GetName() << ": Loaded previous single channel positions from database";
    }
    else {
        LOG(INFO) << GetName() << ": No previous data found";
    }

    // save a copy for comparison at finish stage
    previousAngles = angles;
}



gui::CalibModule_traits::DoFitReturn_t PID_PhiAngle::TheGUI::DoFit(const TH1& hist, unsigned channel)
{
    auto& hist2 = dynamic_cast<const TH2&>(hist);

    h_projection = hist2.ProjectionX("h_projection",channel+1,channel+1);

    func->SetDefaults(h_projection);
    const auto it_fit_param = fitParameters.find(channel);
    if(it_fit_param != fitParameters.end() && !IgnorePreviousFitParameters) {
        VLOG(5) << "Loading previous fit parameters for channel " << channel;
        func->Load(it_fit_param->second);
    }

    for(size_t i=0;i<5;i++)
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

bool PID_PhiAngle::TheGUI::FinishSlice()
{
   h_result = new TGraph(GetNumberOfChannels());

   for(size_t ch=0;ch<GetNumberOfChannels();ch++) {

//       auto angle = angles[ch];

       // always make sure values increase (wrap around)
       /// @todo what to do if slope is negative?
       if(ch!=0)
           while(angles[ch-1] > angles[ch])
               angles[ch] += 360.0;

       h_result->SetPoint(ch, ch, angles[ch]);
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

void PID_PhiAngle::TheGUI::StoreFinishSlice(const interval<TID>& range)
{
    delete h_result;

    TCalibrationData cdata(
                GetName()+"_SingleChannels",
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

