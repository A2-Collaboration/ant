#include "cross_section.h"

#include "base/ParticleType.h"
#include "plot/HistogramFactory.h"
#include "utils/Combinatorics.h"
#include "utils/ParticleTools.h"
#include "base/ParticleTypeTree.h"
#include "slowcontrol/SlowControlVariables.h"
#include "expconfig/ExpConfig.h"

#include "TCanvas.h"

#include <iostream>
#include <memory>


/**
 * Cross section bits
 */

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::utils;

scratch_collicott_CrossSection::Yield_t::Yield_t() {}
scratch_collicott_CrossSection::Scalers_t::Scalers_t() {}


scratch_collicott_CrossSection::scratch_collicott_CrossSection(const HistogramFactory& histFac, OptionsPtr opts) :
    HistFac("CrossSection",histFac)
{
    useSC = opts->Get<bool>("useSC", false);
    if(useSC) {
        slowcontrol::Variables::TaggerScalers->Request();
        slowcontrol::Variables::Trigger->Request();
    }

    Tagger = ExpConfig::Setup::GetDetector<TaggerDetector_t>();
    nTagger = Tagger->GetNChannels();

    flux = HistFac.makeTH1D("flux","Channel","FLUX",BinSettings(nTagger),"flux");

    yield.CreateBranches(HistFac.makeTTree("Yield"));

    scalers.tagger_channel().resize(nTagger);
    scalers.tagger_scalers().resize(nTagger);
    scalers.tagger_eff().resize(nTagger);
    scalers.tagger_deff().resize(nTagger);

    scalers.CreateBranches(HistFac.makeTTree("Scalers"));

}

void scratch_collicott_CrossSection::SetEventType(bool isMC, const string &decay)
{
    event_isMC = isMC;
    event_decay = decay;

    TrackIncidentFlux(); // Track flux if needed
}


void scratch_collicott_CrossSection::AcceptEvent(const LorentzVec& sp, double sp_time, const TTaggerHit& tc, const PromptRandom::Switch &promptrandom)
{
    auto missing = tc.GetPhotonBeam() + LorentzVec(vec3(0,0,0),ParticleTypeDatabase::Proton.Mass()) - sp;

//    h_yield->Fill(sp.Theta()*radtodeg, tc.Channel, promptrandom);
//    missing_mass->Fill(missing.M(), promptrandom);
//    missing_mass_byTC[tc.Channel]->Fill(missing.M(),sp.Theta()*radtodeg, promptrandom);

    yield.isMC      =  event_isMC;
    yield.data      =  event_decay;
    yield.sp_theta  =  sp.Theta()*radtodeg;
    yield.sp_phi    =  sp.Phi()*radtodeg;
    yield.sp_time   =  sp_time;
    yield.sp_im     =  sp.M();

    yield.tc_channel      = tc.Channel;
    yield.tc_photonenergy = tc.PhotonEnergy;
    yield.tc_time         = tc.Time;
    yield.tc_promptrandom = promptrandom.FillWeight();

    yield.missing_mass    = missing.M();

    yield.Tree->Fill();

}

void scratch_collicott_CrossSection::TrackIncidentFlux()
{

    if (!useSC) return;

    // If tagger scalers have not changed, return;
    if(!slowcontrol::Variables::TaggerScalers->HasChanged()) return;

    // Get the new Tagg Scalers
    slowcontrol::Variables::TaggerScalers->GetRates();

    // Fill into tree
    scalers.exp_livetime = slowcontrol::Variables::Trigger->GetExpLivetime();

    for (auto ch = 0u ; ch < nTagger ; ++ch)
    {
        scalers.tagger_channel().at(ch) = ch;
        scalers.tagger_scalers().at(ch) = slowcontrol::Variables::TaggerScalers->GetCounts().at(ch);
        scalers.tagger_eff().at(ch) = Tagger->GetTaggEff(ch).Value;
        scalers.tagger_deff().at(ch) = Tagger->GetTaggEff(ch).Error;

        flux->Fill(ch,scalers.tagger_scalers().at(ch));

    }


    scalers.Tree->Fill();

}
