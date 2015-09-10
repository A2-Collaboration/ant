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
                                            pid_detector->GetNChannels());
}

void PID_Energy::GetGUIs(std::list<std::unique_ptr<gui::Manager_traits> >& guis)
{
    guis.emplace_back(std_ext::make_unique<GUI_Pedestals>(
                          GetName(),
                          Pedestals,
                          calibrationManager,
                          pid_detector
                          ));
}

PID_Energy::ThePhysics::ThePhysics(const string& name,
                                   unsigned nChannels):
    Physics(name)
{
    const BinSettings pid_channels(nChannels);

    h_pedestals = HistFac.makeTH2D(
                      "PID Pedestals",
                      "Raw ADC value",
                      "#",
                      BinSettings(300),
                      pid_channels,
                      "Pedestals");
    h_bananas =
            HistFac.makeTH3D(
                "PID Bananas",
                "CB Energy / MeV",
                "PID Energy / MeV",
                "Channel",
                BinSettings(400,0,800),
                BinSettings(100,0,30),
                pid_channels,
                "Bananas"
                );
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
                h_pedestals->Fill(datum.Value, clusterhit.Channel);
            }
        }
    }

    // bananas
    for(const auto& candidate : event.Reconstructed().Candidates()) {
        if(candidate->Clusters.size() != 2)
            continue;
        if(candidate->Detector() & Detector_t::Type_t::CB &&
           candidate->Detector() & Detector_t::Type_t::PID
           )
        {
            // search for PID cluster
            auto pid_cluster = candidate->FindFirstCluster(Detector_t::Type_t::PID);
            h_bananas->Fill(candidate->ClusterEnergy(),
                            candidate->VetoEnergy(),
                            pid_cluster->CentralElement);
        }
    }
}

void PID_Energy::ThePhysics::Finish()
{
}

void PID_Energy::ThePhysics::ShowResult()
{
    canvas(GetName())
            << drawoption("colz") << h_pedestals
            << drawoption("colz") << h_bananas->Project3D("zy")
            << endc;
}

PID_Energy::GUI_Pedestals::GUI_Pedestals(const string& basename,
                          CalibType& type,
                          const std::shared_ptr<DataManager>& calmgr,
                          const std::shared_ptr<Detector_t>& detector) :
    GUI_CalibType(basename, type, calmgr, detector),
    func(make_shared<gui::FitGausPol0>())
{

}

void PID_Energy::GUI_Pedestals::InitGUI(gui::ManagerWindow_traits* window)
{
    canvas = window->AddCalCanvas();
}

gui::Manager_traits::DoFitReturn_t PID_Energy::GUI_Pedestals::DoFit(TH1* hist, unsigned channel)
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

void PID_Energy::GUI_Pedestals::DisplayFit()
{
    canvas->Show(h_projection, func.get());
}

void PID_Energy::GUI_Pedestals::StoreFit(unsigned channel)
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

bool PID_Energy::GUI_Pedestals::FinishRange()
{
    return false;
}
