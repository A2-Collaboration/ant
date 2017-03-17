#include "TwoPi0_MCSmearing.h"

#include "utils/particle_tools.h"
#include "utils/matcher.h"
#include "utils/Uncertainties.h"
#include "base/Logger.h"
#include "expconfig/ExpConfig.h"
#include "tree/TParticle.h"
#include "tree/TCandidate.h"
#include "tree/TCluster.h"
#include "analysis/utils/uncertainties/FitterSergey.h"
#include "analysis/utils/uncertainties/Interpolated.h"

#include "TH1D.h"
#include "TTree.h"

#include <memory>
#include <cassert>
#include <array>

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace ant::std_ext;
using namespace std;


TwoPi0_MCSmearing::TwoPi0_MCSmearing(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    promptrandom(*ExpConfig::Setup::GetLastFound())
{

    steps = HistFac.makeTH1D("Steps","","#",BinSettings(15),"steps");

    const BinSettings pi0bins(120,80,200);
    const BinSettings thetabins_cb  (35, cos(degree_to_radian(160.0)), cos(degree_to_radian(20.0)));
    const BinSettings thetabins_taps(10, cos(degree_to_radian( 20.0)), cos(degree_to_radian( 0.0)));
    const BinSettings Ebins_theta  (16,0,1600);
    const BinSettings Ebins_element(16,0,1600);

    const auto& cb   = ExpConfig::Setup::GetDetector(Detector_t::Type_t::CB);
    cb_pi0_channel   = HistFac.makeTH2D("CB Pi0",       "m(2#gamma) [MeV]", "Element", pi0bins, BinSettings(cb->GetNChannels()),   "cb_pi0");
    cb_pi0_thetaE    = HistFac.makeTH3D("CB E Theta",   "m(2#gamma) [MeV]", "E_{#gamma} [MeV]", "#cos(#theta)", pi0bins, Ebins_theta, thetabins_cb, "cb_pi0_ETheta");
    cb_pi0_EElement  = HistFac.makeTH3D("CB E element", "m(2#gamma) [MeV]", "E_{#gamma} [MeV]", "Element", pi0bins, Ebins_element, BinSettings(cb->GetNChannels()), "cb_pi0_E_Element");
    cb_channel_correlation = HistFac.makeTH2D("CB Element Correlation",   "Element", "Element", BinSettings(cb->GetNChannels()), BinSettings(cb->GetNChannels()),   "cb_corr");

    const auto& taps = ExpConfig::Setup::GetDetector(Detector_t::Type_t::TAPS);
    taps_pi0_channel  = HistFac.makeTH2D("TAPS Pi0",       "m(2#gamma) [MeV]", "", pi0bins, BinSettings(taps->GetNChannels()), "taps_pi0");
    taps_pi0_thetaE   = HistFac.makeTH3D("TAPS E Theta",   "m(2#gamma) [MeV]", "E_{#gamma} [MeV]", "#cos(#theta)", pi0bins, Ebins_theta, thetabins_taps, "taps_pi0_ETheta");
    taps_pi0_EElement = HistFac.makeTH3D("TAPS E element", "m(2#gamma) [MeV]", "E_{#gamma} [MeV]", "Element", pi0bins, Ebins_element, BinSettings(taps->GetNChannels()), "taps_pi0_E_Element");
    taps_channel_correlation = HistFac.makeTH2D("TAPS Element Correlation",   "Element", "Element", BinSettings(taps->GetNChannels()), BinSettings(taps->GetNChannels()),   "taps_corr");


}

void TwoPi0_MCSmearing::ProcessEvent(const TEvent& event, manager_t&)
{
    const auto& data = event.Reconstructed();

    steps->Fill("Seen",1);

    // cut on energy sum and number of candidates

    if(data.Trigger.CBEnergySum <= 550)
        return;
    steps->Fill("CBESum>550MeV",1);

    const auto& cands = data.Candidates;


    for(auto& c : cands) {
        const auto& cl = c.FindCaloCluster();
        if(cl->HasFlag(TCluster::Flags_t::TouchesHoleCentral))
            return;
    }

    // iterate over tagger hits

    for(const TTaggerHit& taggerhit : data.TaggerHits) {

        steps->Fill("Seen taggerhits",1.0);

        promptrandom.SetTaggerHit(taggerhit.Time);
        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;



    } // Loop taggerhits

}

void TwoPi0_MCSmearing::ShowResult()
{
    canvas(GetName())
            << steps
            << endc;
}

void TwoPi0_MCSmearing::FillIM(const TCandidate& c, double IM)
{
    const auto& cluster = c.FindCaloCluster();

    if(c.Detector & Detector_t::Type_t::CB) {
        cb_pi0_channel->Fill( IM, cluster->CentralElement);
        cb_pi0_thetaE->Fill(  IM, c.CaloEnergy, cos(c.Theta));
        cb_pi0_EElement->Fill(IM, c.CaloEnergy, cluster->CentralElement);
    } else if(c.Detector & Detector_t::Type_t::TAPS) {
        taps_pi0_channel->Fill( IM, cluster->CentralElement);
        taps_pi0_thetaE->Fill(  IM, c.CaloEnergy, cos(c.Theta));
        taps_pi0_EElement->Fill(IM, c.CaloEnergy, cluster->CentralElement);
    }
}

void TwoPi0_MCSmearing::FillCorrelation(const TCandidate& c1, const TCandidate& c2)
{
    if((c1.Detector & Detector_t::Any_t::Calo) != (c2.Detector & Detector_t::Any_t::Calo))
        return;

    const auto& cluster1 = c1.FindCaloCluster();
    const auto& cluster2 = c2.FindCaloCluster();

    if(c1.Detector & Detector_t::Type_t::CB) {
        cb_channel_correlation->Fill(cluster1->CentralElement, cluster2->CentralElement);
        cb_channel_correlation->Fill(cluster2->CentralElement, cluster1->CentralElement);
    } else  if(c1.Detector & Detector_t::Type_t::TAPS) {
        taps_channel_correlation->Fill(cluster1->CentralElement, cluster2->CentralElement);
        taps_channel_correlation->Fill(cluster2->CentralElement, cluster1->CentralElement);
    }
}

AUTO_REGISTER_PHYSICS(TwoPi0_MCSmearing)
