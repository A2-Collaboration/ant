#include "ProtonPi0.h"

#include "utils/uncertainties/Interpolated.h"

#include "base/ParticleTypeTree.h"
#include "base/std_ext/misc.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::physics;

ProtonPi0::ProtonPi0(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    treefitter("treefitter",
               ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0_2g),
               utils::UncertaintyModels::Interpolated::makeAndLoad(),
               true // enable fit z vertex
               )
{
    promptrandom.AddPromptRange({-2.5,2.5});
    promptrandom.AddRandomRange({-50,-5});
    promptrandom.AddRandomRange({  5,50});
    treefitter.SetZVertexSigma(3.0);
    t.CreateBranches(HistFac.makeTTree("t"));

    // set some iteration filter to speed up fitting
    auto pi0 = treefitter.GetTreeNode(ParticleTypeDatabase::Pi0);
    treefitter.SetIterationFilter(
                [pi0] () {
        return ParticleTypeDatabase::Pi0.GetWindow(70).Contains(pi0->Get().LVSum.M());
    });
}

void ProtonPi0::ProcessEvent(const TEvent& event, manager_t&)
{
    const auto& data = event.Reconstructed();
    if(data.Trigger.CBEnergySum<550)
        return;

    if(data.Candidates.size() != 3)
        return;

    // fill some PID info
    t.PID_Ch().clear();
    t.PID_Phi().clear();
    t.PID_E().clear();
    for(const TCluster& cl : data.Clusters) {
        if(cl.DetectorType != Detector_t::Type_t::PID)
            continue;
        t.PID_Ch().push_back(cl.CentralElement);
        t.PID_Phi().push_back(std_ext::radian_to_degree(cl.Position.Phi()));
        t.PID_E().push_back(cl.Energy);
    }


    for(const TTaggerHit& taggerhit : data.TaggerHits) {
        promptrandom.SetTaggerHit(taggerhit.Time);
        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        t.TaggCh = taggerhit.Channel;
        t.TaggE  = taggerhit.PhotonEnergy;
        t.TaggW  = promptrandom.FillWeight();
        t.TaggT  = taggerhit.Time;

        t.FitProb = std_ext::NaN;

        for(auto cand_proton : data.Candidates.get_iter()) {
            auto proton = make_shared<TParticle>(ParticleTypeDatabase::Proton, cand_proton);
            TParticleList photons;
            LorentzVec photon_sum;
            for(auto cand_photon : data.Candidates.get_iter()) {
                if(cand_photon == cand_proton)
                    continue;
                photons.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::Photon, cand_photon));
                photon_sum += *photons.back();
            }

            treefitter.SetEgammaBeam(taggerhit.PhotonEnergy);
            treefitter.SetProton(proton);
            treefitter.SetPhotons(photons);

            APLCON::Result_t fitresult;
            while(treefitter.NextFit(fitresult)) {
                if(fitresult.Status != APLCON::Result_Status_t::Success)
                    continue;
                if(!std_ext::copy_if_greater(t.FitProb, fitresult.Probability))
                    continue;
                auto fitted_proton = treefitter.GetFittedProton();
                t.Proton_Ek = fitted_proton->Ek();
                t.Proton_Phi = std_ext::radian_to_degree(fitted_proton->Phi());
                t.Proton_Theta = std_ext::radian_to_degree(fitted_proton->Theta());
                t.Proton_VetoE = fitted_proton->Candidate->VetoEnergy;
                t.IM_2g = photon_sum.M();
            }
        }

        // require reasonable fit and proton in CB
        if(t.FitProb>0.01 && t.Proton_Theta>25.0) {
            t.Tree->Fill();
        }
    }
}

void ProtonPi0::ShowResult()
{
    canvas(GetName())
            << TTree_drawable(t.Tree, "PID_Phi-Proton_Phi >> (100,-70,70)","")
            << drawoption("colz") << TTree_drawable(t.Tree, "Proton_VetoE:Proton_Ek","")
            << endc;
}




AUTO_REGISTER_PHYSICS(ProtonPi0)
