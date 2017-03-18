#include "gp_pi0[2g]p.h"

#include "base/ParticleType.h"
#include "plot/HistogramFactory.h"
#include "plot/root_draw.h"
#include "utils/combinatorics.h"
#include "utils/particle_tools.h"
#include "base/ParticleTypeTree.h"
#include "utils/uncertainties/FitterSergey.h"
#include "expconfig/ExpConfig.h"
#include "base/Logger.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "base/WrapTFile.h"

#include <iostream>
#include <memory>

/**
 * Single pion photoproduction, decay to 2gamma
 */

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

//ppi0_2gamma::Yield_t::Yield_t() {}

scratch_collicott_ppi0_2gamma::scratch_collicott_ppi0_2gamma(const std::string& name, OptionsPtr opts):
    Physics(name, opts),
    cross_section(HistFac,opts),
    detection_efficiency(HistFac,opts),
    promptrandom(*ExpConfig::Setup::GetLastFound())
{
    steps.CreateBranches(HistFac.makeTTree("analysis_cuts"));
}


void scratch_collicott_ppi0_2gamma::ProcessEvent(const TEvent& event, manager_t&)
{
    triggersimu.ProcessEvent(event);

    // Check the decay string for MC
    // ******************************
    bool signal  = false;
    string decay;
    if(event.Reconstructed().ID.isSet(ant::TID::Flags_t::MC))
            decay = utils::ParticleTools::GetDecayString(event.MCTrue().ParticleTree);
    else    decay = "data" + ExpConfig::Setup::GetLastFound()->GetName();
    if(decay == "(#gamma p) #rightarrow #pi^{0} [ #gamma #gamma ] p ") signal = true;


    steps.AddStep(signal, "All events");

    cross_section.SetEventType(event.Reconstructed().ID.isSet(ant::TID::Flags_t::MC), decay);

    if(event.Reconstructed().ID.isSet(ant::TID::Flags_t::MC))
    detection_efficiency.SetEventType(signal,decay);


    // If MC, track detection efficiency
    // *************************************
    if((event.Reconstructed().ID.isSet(ant::TID::Flags_t::MC)) && (signal)){

        auto pi0 = utils::ParticleTools::FindParticle(ParticleTypeDatabase::Pi0,event.MCTrue().ParticleTree);
        if(!pi0) LOG(ERROR) << "Did not find pi0 in MC True ... strange!!!";

        for (const auto &tc : event.Reconstructed().TaggerHits){

            promptrandom.SetTaggerHit(triggersimu.GetCorrectedTaggerTime(tc));
            detection_efficiency.TrackSignalEvent(*pi0, 0, tc, promptrandom);
        }
    }



    // Get full list of particles, create neutral/charged list
    // *******************************************************
    const auto& candidates = event.Reconstructed().Candidates;
    TCandidatePtrList neutral;
    TCandidatePtrList charged;

    for(const auto& cand : candidates.get_iter()) {
        if(cand->VetoEnergy == 0.0) neutral.emplace_back(cand);
        else charged.emplace_back(cand);
    }



    // ===================================================================
    // Basic event selection
    if (neutral.size() < 2) return; // neutral == 2, 3 is ok
    if (neutral.size() > 3) return; // neutral == 2, 3 is ok
    steps.AddStep(signal, "# Neutral == 2 or 3");
    steps.AddStep(signal, "# Neutral == "+to_string(neutral.size()));

    if (charged.size()  > 1) return; // charged == 0, 1 is ok
    steps.AddStep(signal, "# Charged == 0 or 1");
    steps.AddStep(signal, "# Charged == "+to_string(charged.size()));
    // ===================================================================



    // Create pion candidate
    // *********************
    TParticle Meson;
    double Meson_time = 0;
    for (const auto& photon : neutral)
    {
        Meson+=TParticle(ParticleTypeDatabase::Photon, photon);
        Meson_time+=photon->Time;
    }
    Meson_time = Meson_time / neutral.size(); // Avg the meson time


    // Loop over the Tagger hits
    // *************************
    int i = 0;
    for (const auto &tc : event.Reconstructed().TaggerHits)
    {
        auto missing = tc.GetPhotonBeam() + LorentzVec(vec3(0,0,0),ParticleTypeDatabase::Proton.Mass()) - Meson;

        // ===================================================================
        // Advanced event selection
        if(Meson.M() < 80)  continue;
        if(Meson.M() > 200) continue;
        steps.AddStep(signal, "IM = (80-200)", promptrandom.FillWeight());

        if(charged.size() > 0)
        {
            auto coplanarity  = abs(Meson.Phi() - charged.at(0)->Phi)*radtodeg;
            if (coplanarity < 150) continue;
            if (coplanarity > 210) continue;
            steps.AddStep(signal, "Cop = (150-210)", promptrandom.FillWeight());

        }

        if(missing.M() > 1300) continue;
        steps.AddStep(signal, "MM < 1300", promptrandom.FillWeight());

        // ===================================================================

        promptrandom.SetTaggerHit(triggersimu.GetCorrectedTaggerTime(tc));

        detection_efficiency.AcceptEvent(Meson,Meson_time, tc,promptrandom);
        cross_section.AcceptEvent(Meson,Meson_time,tc,promptrandom);

        i++;
    }


}

void scratch_collicott_ppi0_2gamma::Finish()
{
}

void scratch_collicott_ppi0_2gamma::ShowResult()
{
    ant::canvas(GetName()+": Analysis cuts")
            << TTree_drawable(steps.Tree, "cut.c_str()>>cuts_signal","promptrandom*isSignal")
            << TTree_drawable(steps.Tree, "cut.c_str()>>cuts_background","promptrandom*(!isSignal)")
            << TTree_drawable(steps.Tree, "isSignal","promptrandom")
            << TTree_drawable(steps.Tree, "isSignal","promptrandom","signalcount","0 = bkg, 1 = sig","",BinSettings(2))
            << endc; // actually draws the canvas
}

AUTO_REGISTER_PHYSICS(scratch_collicott_ppi0_2gamma)
