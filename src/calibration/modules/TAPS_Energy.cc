#include "TAPS_Energy.h"

#include "expconfig/detectors/TAPS.h"

#include "analysis/plot/HistogramFactories.h"
#include "analysis/data/Event.h"
#include "analysis/utils/combinatorics.h"

#include "calibration/fitfunctions/FitLandau.h"
#include "calibration/gui/CalCanvas.h"
#include "calibration/fitfunctions/FitGausPol3.h"

#include "tree/TDataRecord.h"

#include "root-addons/cbtaps_display/TH2TAPS.h"
#include "base/Logger.h"

#include <list>
#include <cmath>

using namespace std;
using namespace ant;
using namespace ant::calibration;
using namespace ant::analysis;
using namespace ant::analysis::data;

TAPS_Energy::TAPS_Energy(std::shared_ptr<expconfig::detector::TAPS> taps,
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
           defaultRelativeGain),
    taps_detector(taps)
{

}

TAPS_Energy::ThePhysics::ThePhysics(const string& name, shared_ptr<expconfig::detector::TAPS> taps) :
    Physics(name),
    taps_detector(taps)
{
    const BinSettings taps_channels(taps->GetNChannels());
    const BinSettings energybins(1000);
    const BinSettings timebins(1000,-100,100);

    ggIM = HistFac.makeTH2D("2 neutral IM (TAPS,CB)", "IM [MeV]", "#",
                            energybins, taps_channels, "RelativeGains");
    timing_cuts = HistFac.makeTH2D("Check timing cuts", "IM [MeV]", "#",
                            timebins, taps_channels, "timing_cuts");

    h_pedestals = HistFac.makeTH2D(
                      "TAPS Pedestals",
                      "Raw ADC value",
                      "#",
                      BinSettings(300),
                      taps_channels,
                      "Pedestals");
}

void TAPS_Energy::ThePhysics::ProcessEvent(const Event& event)
{
    const auto& cands = event.Reconstructed.Candidates;

    // pedestals
    for(const Cluster& cluster : event.Reconstructed.AllClusters) {
        if(!(cluster.Detector & Detector_t::Type_t::TAPS))
            continue;
        for(const Cluster::Hit& clusterhit : cluster.Hits) {
            /// \todo check for timing hit?
            /// \todo check for trigger pattern?
            for(const Cluster::Hit::Datum& datum : clusterhit.Data) {
                if(datum.Type != Channel_t::Type_t::Pedestal)
                    continue;
                h_pedestals->Fill(datum.Value, clusterhit.Channel);
            }
        }
    }

    const auto CBTAPS = Detector_t::Type_t::CB | Detector_t::Type_t::TAPS;

    for( auto comb = analysis::utils::makeCombination(cands,2); !comb.Done(); ++comb ) {
        const CandidatePtr& cand1 = comb.at(0);
        const CandidatePtr& cand2 = comb.at(1);

        if(cand1->VetoEnergy==0 && cand2->VetoEnergy==0) {

            //require exactly 1 CB and 1 TAPS
            const auto dets = (cand1->Detector & CBTAPS) ^ (cand2->Detector & CBTAPS);

            if(dets & CBTAPS) {
                const Particle a(ParticleTypeDatabase::Photon,comb.at(0));
                const Particle b(ParticleTypeDatabase::Photon,comb.at(1));
                const TLorentzVector& gg = a + b;


                // Find the one that was in TAPS
                auto cand_taps = cand1->Detector & Detector_t::Type_t::TAPS ? cand1 : cand2;
                auto cl_taps = cand_taps->FindCaloCluster();
                if(cl_taps) {
                    const unsigned ch = cl_taps->CentralElement;
                    const unsigned ring = taps_detector->GetRing(ch);

                    // fill in IM only if ring>4 or if timecut is passed
                    double weight = -1.0;
                    if(ring > 4 || fabs(cand_taps->Time) < 2.5) {
                        weight = 1.0;
                        ggIM->Fill(gg.M(),ch);
                    }
                    timing_cuts->Fill(cand_taps->Time, ch, weight);
                }
            }
        }
    }

}

void TAPS_Energy::ThePhysics::Finish()
{
}

void TAPS_Energy::ThePhysics::ShowResult()
{
    canvas(GetName()) << drawoption("colz") << ggIM
                      << drawoption("colz") << timing_cuts
                      << drawoption("colz") << h_pedestals
                      << endc;
}

unique_ptr<analysis::Physics> TAPS_Energy::GetPhysicsModule()
{
    return std_ext::make_unique<ThePhysics>(GetName(), taps_detector);
}

void TAPS_Energy::GetGUIs(std::list<std::unique_ptr<gui::CalibModule_traits> >& guis)
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




TAPS_Energy::GUI_Gains::GUI_Gains(const string& basename,
                          CalibType& type,
                          const std::shared_ptr<DataManager>& calmgr,
                          const std::shared_ptr<Detector_t>& detector) :
    GUI_CalibType(basename, type, calmgr, detector),
    func(make_shared<gui::FitGausPol3>())
{
}

void TAPS_Energy::GUI_Gains::InitGUI(gui::ManagerWindow_traits* window)
{
    GUI_CalibType::InitGUI(window);

    window->AddNumberEntry("Chi2/NDF limit for autostop", AutoStopOnChi2);

    canvas = window->AddCalCanvas();
    h_peaks = new TH1D("h_peaks","Peak positions",GetNumberOfChannels(),0,GetNumberOfChannels());
    h_peaks->SetXTitle("Channel Number");
    h_peaks->SetYTitle("Pi0 Peak / MeV");
    h_relative = new TH1D("h_relative","Relative change from previous gains",GetNumberOfChannels(),0,GetNumberOfChannels());
    h_relative->SetXTitle("Channel Number");
    h_relative->SetYTitle("Relative change / %");

    h_relative_taps = new TH2TAPS("h_relative_taps",h_relative->GetTitle());
    h_peaks_taps = new TH2TAPS("h_peaks_taps",h_peaks->GetTitle());
}

gui::CalibModule_traits::DoFitReturn_t TAPS_Energy::GUI_Gains::DoFit(TH1* hist, unsigned channel)
{
    if(detector->IsIgnored(channel))
        return DoFitReturn_t::Skip;

    TH2* hist2 = dynamic_cast<TH2*>(hist);

    h_projection = hist2->ProjectionX("h_projection",channel+1,channel+1);

    func->SetDefaults(h_projection);
    func->SetRange(interval<double>(80,250));
    const auto it_fit_param = fitParameters.find(channel);
    if(it_fit_param != fitParameters.end() && !IgnorePreviousFitParameters) {
        VLOG(5) << "Loading previous fit parameters for channel " << channel;
        func->Load(it_fit_param->second);
    }
    else {
        func->FitBackground(h_projection);
    }

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

    if(fit_loop(5))
        return DoFitReturn_t::Next;

    // try with defaults and background fit
    func->SetDefaults(h_projection);
    func->FitBackground(h_projection);

    if(fit_loop(5))
        return DoFitReturn_t::Next;

    // reached maximum retries without good chi2
    LOG(INFO) << "Chi2/dof = " << func->Chi2NDF();
    return DoFitReturn_t::Display;
}

void TAPS_Energy::GUI_Gains::DisplayFit()
{
    canvas->Divide(1,1);
    canvas->Show(h_projection, func.get());
}

void TAPS_Energy::GUI_Gains::StoreFit(unsigned channel)
{
    const double oldValue = previousValues[channel];
    /// \todo obtain convergenceFactor and pi0mass from config or database
    const double convergenceFactor = 1.0;
    const double pi0mass = 135.0;
    const double pi0peak = func->GetPeakPosition();

    // apply convergenceFactor only to the desired procentual change of oldValue,
    // given by (pi0mass/pi0peak - 1)
    const double newValue = oldValue + oldValue * convergenceFactor * (pi0mass/pi0peak - 1);

    calibType.Values[channel] = newValue;

    const double relative_change = 100*(newValue/oldValue-1);

    LOG(INFO) << "Stored Ch=" << channel << ": PeakPosition " << pi0peak
              << " MeV,  gain changed " << oldValue << " -> " << newValue
              << " (" << relative_change << " %)";


    // don't forget the fit parameters
    fitParameters[channel] = func->Save();

    h_peaks->SetBinContent(channel+1, pi0peak);
    h_relative->SetBinContent(channel+1, relative_change);
}

bool TAPS_Energy::GUI_Gains::FinishSlice()
{
    canvas->Clear();
    canvas->Divide(2,2);

    canvas->cd(1);
    h_peaks->SetStats(false);
    h_peaks->Draw("P");
    canvas->cd(2);
    h_peaks_taps->FillElements(*h_peaks);
    h_peaks_taps->Draw("colz");

    canvas->cd(3);
    h_relative->SetStats(false);
    h_relative->Draw("P");
    canvas->cd(4);
    h_relative_taps->FillElements(*h_relative);
    h_relative_taps->Draw("colz");

    return true;
}
