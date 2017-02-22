#include "TestSigmaPlusFitter.h"

#include "base/std_ext/misc.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

#include "utils/uncertainties/Interpolated.h"

TestSigmaPlusFitter::TestSigmaPlusFitter(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    treefitter(name,
               ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::SigmaPlusK0s_6g),
               utils::UncertaintyModels::Interpolated::makeAndLoad()
               )
{

}

void TestSigmaPlusFitter::ProcessEvent(const TEvent& event, manager_t&)
{
    const auto& cands = event.Reconstructed().Candidates;
    if(cands.size() != 7)
        return;

    for(const auto& taggerhit : event.Reconstructed().TaggerHits) {


        for(auto& cand_proton : cands.get_iter()) {

            auto proton = make_shared<TParticle>(ParticleTypeDatabase::Proton, cand_proton);
            TParticleList photons;
            for(auto& cand_photon : cands.get_iter()) {
                if(cand_photon == cand_proton)
                    continue;
                photons.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::Photon, cand_photon));
            }

            treefitter.SetEgammaBeam(taggerhit.PhotonEnergy);
            treefitter.SetProton(proton);
            treefitter.SetPhotons(photons);

            APLCON::Result_t fitresult;
            auto fitprob = std_ext::NaN;
            while(treefitter.NextFit(fitresult)) {
                if(fitresult.Status != APLCON::Result_Status_t::Success)
                    continue;
                if(!std_ext::copy_if_greater(fitprob, fitresult.Probability))
                    continue;
            }
        }

    }


}

void TestSigmaPlusFitter::ShowResult()
{

}

AUTO_REGISTER_PHYSICS(TestSigmaPlusFitter)