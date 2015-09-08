#include "PID_Energy.h"

#include "calibration/gui/CalCanvas.h"
#include "calibration/fitfunctions/FitGausPol0.h"

#include "analysis/plot/HistogramFactories.h"
#include "analysis/data/Event.h"

#include "expconfig/detectors/PID.h"

#include "tree/TDataRecord.h"

#include "base/Logger.h"

#include <list>


using namespace std;
using namespace ant;
using namespace ant::calibration;
using namespace ant::analysis;
using namespace ant::analysis::data;

PID_Energy::PID_Energy(
        std::shared_ptr<expconfig::detector::PID> pid,
        std::shared_ptr<DataManager> calmgr,
        Calibration::Converter::ptr_t converter,
        double defaultPedestal,
        double defaultGain,
        double defaultThreshold,
        double defaultRelativeGain
        ) :
    Energy(pid->Type,
           calmgr,
           converter,
           defaultPedestal,
           defaultGain,
           defaultThreshold,
           defaultRelativeGain),
    pid_detector(pid)
{

}


PID_Energy::~PID_Energy()
{

}

unique_ptr<analysis::Physics> PID_Energy::GetPhysicsModule()
{
    return std_ext::make_unique<ThePhysics>(GetName(),
                                            Pedestals.Name,
                                            pid_detector->GetNChannels());
}

void PID_Energy::GetGUIs(std::list<std::unique_ptr<gui::Manager_traits> >& guis)
{
    guis.emplace_back(std_ext::make_unique<TheGUI>(
                          GetName(),
                          Pedestals,
                          calibrationManager,
                          pid_detector
                          ));
}

PID_Energy::ThePhysics::ThePhysics(const string& name,
                                   const string& hist_name,
                                   unsigned nChannels):
    Physics(name)
{
    const BinSettings pid_channels(nChannels);

    pedestals = HistFac.makeTH2D("PID Pedestals", "Raw ADC value", "#", BinSettings(300), pid_channels, hist_name);
}

void PID_Energy::ThePhysics::ProcessEvent(const data::Event& event)
{

    // pedestals, best determined from insane clusters (that is no timing hit seen)
    for(const Cluster& cluster : event.Reconstructed().InsaneClusters()) {
        if(!(cluster.Detector & Detector_t::Type_t::PID))
            continue;
        for(const Cluster::Hit& clusterhit : cluster.Hits) {
            for(const Cluster::Hit::Datum& datum : clusterhit.Data) {
                if(datum.Type != Channel_t::Type_t::Pedestal)
                    continue;
                pedestals->Fill(datum.Value, clusterhit.Channel);
            }
        }
    }
}

void PID_Energy::ThePhysics::Finish()
{
}

void PID_Energy::ThePhysics::ShowResult()
{
    canvas(GetName()) << drawoption("colz") << pedestals << endc;
}

PID_Energy::TheGUI::TheGUI(const string& basename,
                          CalibType& type,
                          const std::shared_ptr<DataManager>& calmgr,
                          const std::shared_ptr<Detector_t>& detector) :
    GUI_CalibType(basename, type, calmgr, detector),
    func(make_shared<gui::FitGausPol0>())
{

}

void PID_Energy::TheGUI::InitGUI(gui::ManagerWindow_traits* window)
{
    canvas = window->AddCalCanvas();
}

gui::Manager_traits::DoFitReturn_t PID_Energy::TheGUI::DoFit(TH1* hist, unsigned channel)
{
    if(detector->IsIgnored(channel))
        return DoFitReturn_t::Skip;

    TH2* hist2 = dynamic_cast<TH2*>(hist);

    h_projection = hist2->ProjectionX("",channel+1,channel+1);

    func->SetDefaults(h_projection);
    const auto it_fit_param = fitParameters.find(channel);
    if(it_fit_param != fitParameters.end()) {
        VLOG(5) << "Loading previous fit parameters for channel " << channel;
        func->Load(it_fit_param->second);
    }

    func->Fit(h_projection);

    /// \todo implement automatic stop if fit failed?

    // goto next channel
    return DoFitReturn_t::Next;
}

void PID_Energy::TheGUI::DisplayFit()
{
    canvas->Show(h_projection, func.get());
}

void PID_Energy::TheGUI::StoreFit(unsigned channel)
{

    const double oldValue = previousValues[channel];
    const double newValue = func->GetPeakPosition();

    calibType.Values[channel] = newValue;

    const double relative_change = 100*(newValue/oldValue-1);

    LOG(INFO) << "Stored Ch=" << channel << ":  "
              <<" Pedestal changed " << oldValue << " -> " << newValue
              << " (" << relative_change << " %)";


    // don't forget the fit parameters
    fitParameters[channel] = func->Save();
}

bool PID_Energy::TheGUI::FinishRange()
{
    return false;
}
