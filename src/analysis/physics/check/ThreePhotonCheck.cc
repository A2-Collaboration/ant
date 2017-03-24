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
    h_photonsIM = HistFac.makeTH1D("3#gamma IM","MeV","",BinSettings(400),"h_PhotonIM");
}

struct MaxTracker {
    double v;
    MaxTracker(const double& start = std_ext::inf) : v(start) {}
    bool Track(const double& value) { if(value > v) {v = value; return true;} else return false; }
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

    const auto nTAPS = count_if(cands.begin(), cands.end(),
                                [] (const TCandidate& c) { return c.Detector & Detector_t::Type_t::TAPS;});

    if(nTAPS > 1)
        return;

    h_Steps->Fill("nTAPS==0,1",1.0);


    for(const TTaggerHit& taggerhit : data.TaggerHits) {

        promptrandom.SetTaggerTime(triggersimu.GetCorrectedTaggerTime(taggerhit));
        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;
        h_Steps->Fill("TagHits",1.0);
        if(promptrandom.State() == PromptRandom::Case::Prompt)
            h_Steps->Fill("TagHits prompt",1.0);

        MaxTracker best_tacker(-1.0);
        LorentzVec best_photon_sum;

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

            const auto res = fitter.DoFit(taggerhit.PhotonEnergy, proton, photons);

            if(best_tacker.Track(res.Probability)) {
                best_photon_sum = photon_sum;
            }
        }

        if(isfinite(best_tacker()) && best_tacker() > 0.05) {
            h_Steps->Fill("Fit Ok",1.0);
            h_photonsIM->Fill(best_photon_sum.M(), promptrandom.FillWeight());
        }

    }
}

void ThreePhotonCheck::ShowResult()
{
    canvas(GetName())
            << h_Steps
            << h_photonsIM
            << endc;
}




AUTO_REGISTER_PHYSICS(ThreePhotonCheck)
