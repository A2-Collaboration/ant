
#include "CB_SourceCalib.h"

#include "calibration/gui/CalCanvas.h"
#include "calibration/fitfunctions/FitGausexpo.h"

#include "expconfig/detectors/CB.h"


#include "root-addons/cbtaps_display/TH2CB.h"
#include "base/Logger.h"

using namespace std;
using namespace ant;
using namespace ant::calibration;

CB_SourceCalib::CB_SourceCalib(std::shared_ptr<expconfig::detector::CB> cb,
                               std::shared_ptr<DataManager> calmgr,
                               Calibration::Converter::ptr_t converter,
                               double defaultPedestal,
                               double defaultGain,
                               double defaultThreshold,
                               double defaultRelativeGain):
    Energy(cb->Type,
           calmgr,
           converter,
           defaultPedestal,
           defaultGain,
           defaultThreshold,
           defaultRelativeGain),
    cb_detector(cb)
{

}

void CB_SourceCalib::GetGUIs(list<unique_ptr<gui::CalibModule_traits> >& guis)
{
    guis.emplace_back(std_ext::make_unique<GUI_Gains>(
                          GetName(),
                          RelativeGains,
                          calibrationManager,
                          cb_detector
                          ));
}

CB_SourceCalib::GUI_Gains::GUI_Gains(const string& basename,
                                     CalibType& type,
                                     const std::shared_ptr<DataManager>& calmgr,
                                     const std::shared_ptr<Detector_t>& detector) :
               GUI_CalibType(basename, type, calmgr, detector),
               func(make_shared<gui::FitGausexpo>())
{
}


void CB_SourceCalib::GUI_Gains::InitGUI(gui::ManagerWindow_traits* window)
{
    GUI_CalibType::InitGUI(window);

    window->AddNumberEntry("Chi2/NDF limit for autostop", AutoStopOnChi2);

    canvas = window->AddCalCanvas();
    h_peaks = new TH1D("h_peaks","Peak positions",GetNumberOfChannels(),0,GetNumberOfChannels());
    h_peaks-> SetXTitle("Channel Number");
    h_peaks-> SetYTitle("4.4Mev Peak");
}

gui::CalibModule_traits::DoFitReturn_t CB_SourceCalib::GUI_Gains::DoFit(TH1 *hist, unsigned channel)
{
    if(detector->IsIgnored(channel))
        return DoFitReturn_t::Skip;

    TH2* hist2 = dynamic_cast<TH2*>(hist);

    h_projection = hist2->ProjectionX("h_projection",channel+1,channel+1);

    if(h_projection->GetEntries()==0)
        return DoFitReturn_t::Display;

    func->SetDefaults(h_projection);

    const auto it_fit_param =fitParameters.find(channel);

    auto fit_loop = [this] (size_t retries) {
        do {
            func->Fit(h_projection);
            VLOG(5) << "Chi2/dof = " << func->Chi2NDF();
            if(func->Chi2NDF() < AutoStopOnChi2) {
                return true;
            }
            retries--;
        }
        while(retries>0);
        return false;
    };


}
void CB_SourceCalib::GUI_Gains::DisplayFit()
{
    canvas->Divide(1,1);
    canvas->Show(h_projection, func.get());
}

void CB_SourceCalib::GUI_Gains::StoreFit(unsigned channel)
{
    const double oldValue = previousValues[channel];


}

shared_ptr<TH1> CB_SourceCalib::GUI_Gains::GetHistogram(const WrapTFile &file) const
{
    return file.GetSharedHist<TH1>("CB_SourceCalib/HitsADC_Cluster");
}


bool CB_SourceCalib::GUI_Gains::FinishSlice()
{
    canvas->Clear();
    canvas->Divide(2,2);

    canvas->cd(1);
    h_peaks->SetStats(false);
    h_peaks->Draw("P");

    return true;
}

string CB_SourceCalib::GUI_Gains::GetName() const
{
    return "bla";
}










