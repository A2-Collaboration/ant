#include "CB_Energy.h"

#include "calibration/gui/CalCanvas.h"

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
    p(parent)
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
    return false;
}

void CB_Energy::TheGUI::DisplayFit()
{
    LOG(INFO) << "Displaying Fit";
}

void CB_Energy::TheGUI::StoreFit(unsigned channel)
{
    LOG(INFO) << "Storing data for channel " << channel;
}

bool CB_Energy::TheGUI::FinishRange()
{
    LOG(INFO) << "FinishRange";
    return false;
}
