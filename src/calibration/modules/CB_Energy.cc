#include "CB_Energy.h"

#include "calibration/gui/CalCanvas.h"
#include "calibration/gui/FitGausPol3.h"

#include "analysis/plot/HistogramFactories.h"
#include "analysis/data/Event.h"
#include "analysis/utils/combinatorics.h"

#include "expconfig/detectors/CB.h"

#include "tree/TDataRecord.h"

#include "base/Logger.h"

#include <list>


using namespace std;
using namespace ant;
using namespace ant::calibration;

CB_Energy::CB_Energy(std::shared_ptr<expconfig::detector::CB> cb,
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

unique_ptr<Physics> CB_Energy::GetPhysicsModule()
{
    return std_ext::make_unique<ThePhysics>(GetName(),
                                            GUI_CalibType::ConstructName(GetName(), Gains.Name),
                                            cb_detector->GetNChannels());
}

void CB_Energy::GetGUIs(list<unique_ptr<gui::Manager_traits> >& guis) {
    guis.emplace_back(std_ext::make_unique<TheGUI>(GetName(), Gains, this));
}


CB_Energy::ThePhysics::ThePhysics(const string& name,
                                  const string& hist_name,
                                  unsigned nChannels):
    Physics(name)
{
    const BinSettings cb_channels(nChannels);
    const BinSettings energybins(1000);

    ggIM = HistFac.makeTH2D("2 neutral IM (CB,CB)", "IM [MeV]", "#",
                            energybins, cb_channels, hist_name);
}

void CB_Energy::ThePhysics::ProcessEvent(const Event& event)
{
    const auto& cands = event.Reconstructed().Candidates();

    for( auto comb = makeCombination(cands,2); !comb.Done(); ++comb ) {
        const CandidatePtr& p1 = comb.at(0);
        const CandidatePtr& p2 = comb.at(1);

        if(p1->VetoEnergy()==0 && p2->VetoEnergy()==0
           && (p1->Detector() & Detector_t::Type_t::CB)
           && (p2->Detector() & Detector_t::Type_t::CB)) {
            const Particle a(ParticleTypeDatabase::Photon,comb.at(0));
            const Particle b(ParticleTypeDatabase::Photon,comb.at(1));
            const TLorentzVector gg = a + b;

            auto cl1 = p1->FindCaloCluster();
            if(cl1)
                ggIM->Fill(gg.M(),cl1->CentralElement);

            auto cl2 = p2->FindCaloCluster();
            if(cl2)
                ggIM->Fill(gg.M(),cl2->CentralElement);
        }
    }
}

void CB_Energy::ThePhysics::Finish()
{

}

void CB_Energy::ThePhysics::ShowResult()
{
    canvas(GetName()) << drawoption("colz") << ggIM << endc;
}

CB_Energy::TheGUI::TheGUI(const string& basename, CalibType& type, CB_Energy* parent) :
    GUI_CalibType(basename, type, p->calibrationManager),
    p(parent),
    func(make_shared<gui::FitGausPol3>())
{

}

unsigned CB_Energy::TheGUI::GetNumberOfChannels() const
{
    return p->cb_detector->GetNChannels();
}

void CB_Energy::TheGUI::InitGUI()
{
    c_fit = new gui::CalCanvas(GetName()+": Fit");
    c_overview = new gui::CalCanvas(GetName()+": Overview");
}

list<gui::CalCanvas*> CB_Energy::TheGUI::GetCanvases() const
{
    return {c_fit, c_overview};
}

bool CB_Energy::TheGUI::DoFit(TH1* hist, unsigned channel)
{
    TH2* hist2 = dynamic_cast<TH2*>(hist);

    if(hist2) {
        projection = hist2->ProjectionX("",channel,channel+1);

        func->SetDefaults(projection);
        const auto it_fit_param = fitParameters.find(channel);
        if(it_fit_param != fitParameters.end()) {
            VLOG(5) << "Loading previous fit parameters for channel " << channel;
            func->Load(it_fit_param->second);
        }

        func->Fit(projection);

        if(channel==0) {
            return true;
        }

    } else {
        LOG(WARNING) << "Supplied Hist is not 2D";
    }

    // do not show something
    return false;
}

void CB_Energy::TheGUI::DisplayFit()
{
    c_fit->Show(projection, func.get());
}

void CB_Energy::TheGUI::StoreFit(unsigned channel)
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

    LOG(INFO) << "Stored Ch=" << channel << ": PeakPosition " << pi0peak
              << " MeV,  gain changed " << oldValue << " -> " << newValue
              << " (" << 100*newValue/oldValue << " %)";


}

bool CB_Energy::TheGUI::FinishRange()
{
//    c_overview->Clear();
//    c_overview->cd();
//    TH1D* overview = new TH1D("Overview", "");
//    c_overview->Draw();

    LOG(INFO) << "FinishRange";
    return false;
}

void CB_Energy::TheGUI::StoreFinishRange(const interval<TID>& range)
{
    // c_overview->Close();

    GUI_CalibType::StoreFinishRange(range);
}
