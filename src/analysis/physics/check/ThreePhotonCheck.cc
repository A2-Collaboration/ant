#include "ThreePhotonCheck.h"

#include "utils/uncertainties/Interpolated.h"

#include "base/ParticleTypeTree.h"
#include "base/std_ext/misc.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::physics;

ThreePhotonCheck::ThreePhotonCheck(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    fitter(utils::UncertaintyModels::Interpolated::makeAndLoad(),
           true // enable fit z vertex
           )
{
    promptrandom.AddPromptRange({-3,3});
    promptrandom.AddRandomRange({-20,-10});
    promptrandom.AddRandomRange({ 10,20});
    fitter.SetZVertexSigma(3.0);

    h_Steps     = HistFac.makeTH1D("Steps","","",BinSettings(10),"h_Steps");
    h_photonsIM_raw = HistFac.makeTH1D("3#gamma IM raw","MeV","",{200,700.0,900.0},"h_PhotonIM_raw");
    h_photonsIM_fitted = HistFac.makeTH1D("3#gamma IM fitted","MeV","",{200,700.0,900.0},"h_PhotonIM_fitted");
    h_mm = HistFac.makeTH1D("mm","MeV","",{400,800,1200.0},"h_mm");

}

struct MinTracker {
    double v;
    MinTracker(const double& start = std_ext::inf) : v(start) {}
    bool Track(const double& value) { if(value < v) {v = value; return true;} else return false; }
    double operator()() const { return v; }
};

void ThreePhotonCheck::ProcessEvent(const TEvent& event, manager_t&)
{
    triggersimu.ProcessEvent(event);

    h_Steps->Fill("Seen",1.0);

    const auto& data = event.Reconstructed();
    const auto& cands = data.Candidates;

    if(!triggersimu.HasTriggered())
        return;
    h_Steps->Fill("Triggered",1.0);


    if(cands.size() != 4)
        return;
    h_Steps->Fill("nCands==4",1.0);

     for(const TTaggerHit& taggerhit : data.TaggerHits) {

        promptrandom.SetTaggerTime(triggersimu.GetCorrectedTaggerTime(taggerhit));
        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;
        h_Steps->Fill("TagHits",1.0);
        if(promptrandom.State() == PromptRandom::Case::Prompt)
            h_Steps->Fill("TagHits prompt",1.0);

        MinTracker best_tacker;
        LorentzVec best_photon_sum_raw;
        LorentzVec best_photon_sum_fitted;

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

            const LorentzVec beam_target = taggerhit.GetPhotonBeam() + LorentzVec({},ParticleTypeDatabase::Proton.Mass());
            const auto mm = beam_target - photon_sum - TParticle(ParticleTypeDatabase::Proton, cand_proton);

            if(best_tacker.Track(std_ext::sqr(ParticleTypeDatabase::Proton.Mass() - mm.M()))) {
                best_photon_sum_raw = photon_sum;

                const auto res = fitter.DoFit(taggerhit.PhotonEnergy, proton, photons);

                if(res.Status == APLCON::Result_Status_t::Success) {
                    const auto fitted_photons = fitter.GetFittedPhotons();
                    best_photon_sum_fitted = LVSum(fitted_photons.begin(), fitted_photons.end());
                }
            }
        }

        if(isfinite(best_tacker())) {
            h_photonsIM_raw->Fill(best_photon_sum_raw.M(), promptrandom.FillWeight());
            h_photonsIM_fitted->Fill(best_photon_sum_fitted.M(), promptrandom.FillWeight());
        }

    }
}

void ThreePhotonCheck::ShowResult()
{
    canvas(GetName())
            << h_Steps
            << h_photonsIM_raw
            << h_photonsIM_fitted
            << h_mm
            << endc;
}




AUTO_REGISTER_PHYSICS(ThreePhotonCheck)
