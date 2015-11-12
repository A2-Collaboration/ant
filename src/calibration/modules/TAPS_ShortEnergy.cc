#include "TAPS_ShortEnergy.h"
#include "TF1.h"
#include "expconfig/detectors/TAPS.h"
#include "calibration/fitfunctions/FitLandau.h"
#include "analysis/plot/HistogramFactories.h"
#include "analysis/data/Event.h"
#include "analysis/utils/combinatorics.h"

#include "tree/TDataRecord.h"

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
                      "relative");
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
                        double central_e = 0.0;
                        for(const Cluster::Hit::Datum& datum : clusterhit.Data) {

                            if(datum.Type == Channel_t::Type_t::Integral)
                                central_e = datum.Value;

                            if(datum.Type == Channel_t::Type_t::IntegralShort)
                                h_rel_gamma->Fill(datum.Value / central_e, clusterhit.Channel);

                        }
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



void TAPS_ShortEnergy::GetGUIs(std::list<std::unique_ptr<gui::Manager_traits> >& guis)
{
    guis.emplace_back(std_ext::make_unique<SGGUI_Pedestals>(
                          GetName(),
                          Pedestals,
                          calibrationManager,
                          taps_detector
                          ));
}

// sigma = 3
// A  =
struct myLandau : public gui::FitLandau {
    using FitLandau::FitLandau;

    virtual void SetDefaults(TH1* hist) override {
        func->SetParameter(2, 3.0);
        SetRange({0, 200});

        if(hist) {
            func->SetParameter(0, hist->GetMaximum()/3.0);
            const double max_pos = hist->GetXaxis()->GetBinCenter(hist->GetMaximumBin());
            func->SetParameter(1,max_pos);

        } else {
            func->SetParameter(0,1000);
            func->SetParameter(1,100);
        }
    }

    virtual ~myLandau();
};

myLandau::~myLandau() {}

TAPS_ShortEnergy::SGGUI_Pedestals::SGGUI_Pedestals(
        const string& basename,
        Energy::CalibType& type,
        const std::shared_ptr<DataManager>& calmgr,
        const std::shared_ptr<Detector_t>& detector):
    Energy::GUI_Pedestals::GUI_Pedestals(basename, type, calmgr, detector)
{
    func = make_shared<myLandau>();
    //TODO: subclass fitLandau with specific settings
}
