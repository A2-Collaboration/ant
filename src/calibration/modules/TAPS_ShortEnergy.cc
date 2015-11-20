#include "TAPS_ShortEnergy.h"
#include "TF1.h"
#include "expconfig/detectors/TAPS.h"
#include "calibration/fitfunctions/FitLandau.h"
#include "calibration/fitfunctions/FitGaus.h"
#include "analysis/plot/HistogramFactories.h"
#include "analysis/data/Event.h"
#include "analysis/utils/combinatorics.h"
#include "calibration/gui/CalCanvas.h"
#include "tree/TDataRecord.h"
#include "base/Logger.h"
#include "base/cbtaps_display/TH2TAPS.h"

#include <list>
#include <cmath>

using namespace std;
using namespace ant;
using namespace ant::calibration;
using namespace ant::analysis;
using namespace ant::analysis::data;

TAPS_ShortEnergy::TAPS_ShortEnergy(std::shared_ptr<expconfig::detector::TAPS> taps,
        std::shared_ptr<DataManager> calmgr,
        Calibration::Converter::ptr_t converter,
        double defaultPedestal,
        double defaultGain,
        double defaultThreshold,
        double defaultRelativeGain) :
    Energy(Detector_t::Type_t::TAPS,
           calmgr,
           converter,
           defaultPedestal,
           defaultGain,
           defaultThreshold,
           defaultRelativeGain,
           Channel_t::Type_t::IntegralShort),
    taps_detector(taps)
{

}

TAPS_ShortEnergy::ThePhysics::ThePhysics(const string& name, shared_ptr<expconfig::detector::TAPS> taps) :
    Physics(name),
    taps_detector(taps)
{
    const BinSettings taps_channels(taps->GetNChannels());

    h_pedestals = HistFac.makeTH2D(
                      "TAPS ShortGate Pedestals",
                      "Raw ADC value",
                      "#",
                      BinSettings(300),
                      taps_channels,
                      "Pedestals");

    h_rel_gamma = HistFac.makeTH2D(
                      "TAPS E_{S} / E_{L}",
                      "rel",
                      "#",
                      BinSettings(400,-10,10),
                      taps_channels,
                      "RelativeGains");
}

void TAPS_ShortEnergy::ThePhysics::ProcessEvent(const Event& event)
{
    // pedestals
    for(const Cluster& cluster : event.Reconstructed().AllClusters()) {
        if(!(cluster.Detector & Detector_t::Type_t::TAPS))
            continue;
        for(const Cluster::Hit& clusterhit : cluster.Hits) {
            /// \todo check for timing hit?
            /// \todo check for trigger pattern?
            for(const Cluster::Hit::Datum& datum : clusterhit.Data) {

                if(datum.Type != Channel_t::Type_t::PedestalShort)
                    continue;
                h_pedestals->Fill(datum.Value, clusterhit.Channel);

            }
        }
    }

    for(const auto& c : event.Reconstructed().Candidates()) {
        if(c->VetoEnergy() < 0.5) {
            const auto& cluster = c->FindCaloCluster();

            if(cluster)
                for(const Cluster::Hit& clusterhit : cluster->Hits) {

                    if(clusterhit.Channel == cluster->CentralElement) {
                        double central_e = numeric_limits<double>::quiet_NaN();
                        double short_e = numeric_limits<double>::quiet_NaN();
                        for(const Cluster::Hit::Datum& datum : clusterhit.Data) {

                            if(datum.Type == Channel_t::Type_t::Integral)
                                central_e = datum.Value;

                            if(datum.Type == Channel_t::Type_t::IntegralShort)
                                short_e = datum.Value;
                        }

                        h_rel_gamma->Fill(short_e / central_e, clusterhit.Channel);

                    }
                }
        }
    }
}

void TAPS_ShortEnergy::ThePhysics::Finish()
{
}

void TAPS_ShortEnergy::ThePhysics::ShowResult()
{
    canvas(GetName()) << drawoption("colz") << h_pedestals << h_rel_gamma
                      << endc;
}

unique_ptr<analysis::Physics> TAPS_ShortEnergy::GetPhysicsModule()
{
    return std_ext::make_unique<ThePhysics>(GetName(), taps_detector);
}

void TAPS_ShortEnergy::GetGUIs(std::list<std::unique_ptr<gui::CalibModule_traits> >& guis)
{
    guis.emplace_back(std_ext::make_unique<GUI_Pedestals>(
                          GetName(),
                          Pedestals,
                          calibrationManager,
                          taps_detector,
                          make_shared<gui::FitLandau>()
                          ));
    guis.emplace_back(std_ext::make_unique<GUI_Gains>(
                          GetName(),
                          RelativeGains,
                          calibrationManager,
                          taps_detector
                          ));
}


TAPS_ShortEnergy::GUI_Gains::GUI_Gains(const string& basename,
                          CalibType& type,
                          const std::shared_ptr<DataManager>& calmgr,
                          const std::shared_ptr<expconfig::detector::TAPS>& taps) :
    GUI_CalibType(basename, type, calmgr, taps),
    func(make_shared<gui::FitGaus>()),
    taps_detector(taps)
{

}

TAPS_ShortEnergy::GUI_Gains::~GUI_Gains()
{

}

void TAPS_ShortEnergy::GUI_Gains::InitGUI(gui::ManagerWindow_traits* window)
{
    canvas = window->AddCalCanvas();
    h_peaks = new TH1D("h_peaks","Peak positions",GetNumberOfChannels(),0,GetNumberOfChannels());
    h_peaks->SetXTitle("Channel Number");
    h_peaks->SetYTitle("Peak");
    h_relative = new TH1D("h_relative","Relative change from previous gains",GetNumberOfChannels(),0,GetNumberOfChannels());
    h_relative->SetXTitle("Channel Number");
    h_relative->SetYTitle("Relative change / %");
}

gui::CalibModule_traits::DoFitReturn_t TAPS_ShortEnergy::GUI_Gains::DoFit(TH1* hist, unsigned channel,
                                                               const CalibModule_traits::DoFitOptions_t& options)
{
    if(detector->IsIgnored(channel) || taps_detector->IsPbWO4(channel))
        return DoFitReturn_t::Skip;

    TH2* hist2 = dynamic_cast<TH2*>(hist);

    h_projection = hist2->ProjectionX("",channel+1,channel+1);

    func->SetDefaults(h_projection);
    func->SetRange(interval<double>(-1,3));
    const auto it_fit_param = fitParameters.find(channel);
    if(it_fit_param != fitParameters.end()
       && !options.IgnorePreviousFitParameters) {
        VLOG(5) << "Loading previous fit parameters for channel " << channel;
        func->Load(it_fit_param->second);
    }

    for(size_t i=0;i<5;i++)
        func->Fit(h_projection);

    /// \todo implement automatic stop if fit failed?

    // goto next channel
    return DoFitReturn_t::Next;
}

void TAPS_ShortEnergy::GUI_Gains::DisplayFit()
{
    canvas->Divide(1,1);
    canvas->Show(h_projection, func.get());
}

void TAPS_ShortEnergy::GUI_Gains::StoreFit(unsigned channel)
{
    const double oldValue = previousValues[channel];

    const double convergenceFactor = 1.0;
    const double target_pos = 1.0;
    const double peak = func->GetPeakPosition();

    // apply convergenceFactor only to the desired procentual change of oldValue,
    // given by (pi0mass/pi0peak - 1)
    const double newValue = oldValue + oldValue * convergenceFactor * (target_pos/peak - 1);

    calibType.Values[channel] = newValue;

    const double relative_change = 100*(newValue/oldValue-1);

    LOG(INFO) << "Stored Ch=" << channel << ": PeakPosition " << peak
              << ", gain changed " << oldValue << " -> " << newValue
              << " (" << relative_change << " %)";


    // don't forget the fit parameters
    fitParameters[channel] = func->Save();

    h_peaks->SetBinContent(channel+1, peak);
    h_relative->SetBinContent(channel+1, relative_change);
}

bool TAPS_ShortEnergy::GUI_Gains::FinishSlice()
{
    canvas->Clear();
    canvas->Divide(2,2);

    canvas->cd(1);
    h_peaks->SetStats(false);
    h_peaks->Draw("P");
    canvas->cd(2);
    TH2TAPS* h_peaks_taps = new TH2TAPS("h_peaks_taps",h_peaks->GetTitle());
    h_peaks_taps->FillElements(*h_peaks);
    h_peaks_taps->Draw("colz");

    canvas->cd(3);
    h_relative->SetStats(false);
    h_relative->Draw("P");
    canvas->cd(4);
    TH2TAPS* h_relative_taps = new TH2TAPS("h_relative_cb",h_relative->GetTitle());
    h_relative_taps->FillElements(*h_relative);
    h_relative_taps->Draw("colz");

    return true;
}


TAPS_ShortEnergy::GUI_Pedestals::GUI_Pedestals(const string& basename,
                                               Energy::CalibType& type,
                                               const std::shared_ptr<DataManager>& calmgr,
                                               const std::shared_ptr<expconfig::detector::TAPS>& taps,
                                               std::shared_ptr<gui::PeakingFitFunction> fitfunction) :
    Energy::GUI_Pedestals(basename, type, calmgr, taps, fitfunction),
    taps_detector(taps)
{
}

gui::CalibModule_traits::DoFitReturn_t TAPS_ShortEnergy::GUI_Pedestals::DoFit(TH1* hist, unsigned channel, const gui::CalibModule_traits::DoFitOptions_t& options)
{
    if(taps_detector->IsPbWO4(channel))
        return DoFitReturn_t::Skip;
    return Energy::GUI_Pedestals::DoFit(hist, channel, options);
}
