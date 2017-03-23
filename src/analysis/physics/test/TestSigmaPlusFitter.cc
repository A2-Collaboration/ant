#include "TestSigmaPlusFitter.h"

#include "utils/ParticleTools.h"

#include "base/Logger.h"
#include "base/std_ext/misc.h"


using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

#include "utils/uncertainties/Interpolated.h"

TestSigmaPlusFitter::TestSigmaPlusFitter(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    treefitter(ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::SigmaPlusK0s_6g),
               utils::UncertaintyModels::Interpolated::makeAndLoad()
               ),
    treefitter_SigmaPlus(treefitter.GetTreeNode(ParticleTypeDatabase::SigmaPlus)),
    treefitter_K0s(treefitter.GetTreeNode(ParticleTypeDatabase::K0s))
{
    t.CreateBranches(HistFac.makeTTree("t"));

    // be lazy and catch complete class...
    treefitter.SetIterationFilter([this] () {
        const auto sigmaPlus_cut = ParticleTypeDatabase::SigmaPlus.GetWindow(200);
        const auto K0s_cut = ParticleTypeDatabase::K0s.GetWindow(100);
        auto ok = sigmaPlus_cut.Contains(treefitter_SigmaPlus->Get().LVSum.M()) &&
                  K0s_cut.Contains(treefitter_K0s->Get().LVSum.M());
        return ok;
    });
}

void TestSigmaPlusFitter::ProcessEvent(const TEvent& event, manager_t&)
{

    const auto& cands = event.Reconstructed().Candidates;

    if(cands.size() != 7)
        return;

    for(const auto& taggerhit : event.Reconstructed().TaggerHits) {

        auto fitprob = std_ext::NaN;
        LorentzVec best_SigmaPlus;
        LorentzVec best_K0s;

        for(auto& cand_proton : cands.get_iter()) {

            auto proton = make_shared<TParticle>(ParticleTypeDatabase::Proton, cand_proton);
            TParticleList photons;
            for(auto& cand_photon : cands.get_iter()) {
                if(cand_photon == cand_proton)
                    continue;
                photons.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::Photon, cand_photon));
            }

            treefitter.PrepareFits(taggerhit.PhotonEnergy, proton, photons);

            APLCON::Result_t fitresult;
            while(treefitter.NextFit(fitresult)) {
                if(fitresult.Status != APLCON::Result_Status_t::Success)
                    continue;
                if(!std_ext::copy_if_greater(fitprob, fitresult.Probability))
                    continue;

                best_SigmaPlus = treefitter_SigmaPlus->Get().LVSum;
                best_K0s = treefitter_K0s->Get().LVSum;
            }
        }

        if(fitprob > 0.01 && event.MCTrue().ParticleTree) {
            auto mctree = event.MCTrue().ParticleTree;
            auto true_SigmaPlus = utils::ParticleTools::FindParticle(ParticleTypeDatabase::SigmaPlus, mctree);
            auto true_K0s = utils::ParticleTools::FindParticle(ParticleTypeDatabase::K0s, mctree);

            t.SigmaPlus_DeltaAngle = std_ext::radian_to_degree(true_SigmaPlus->Angle(best_SigmaPlus));
            t.SigmaPlus_DeltaE     = (true_SigmaPlus->E - best_SigmaPlus.E)/true_SigmaPlus->E;

            t.K0s_DeltaAngle = std_ext::radian_to_degree(true_K0s->Angle(best_K0s));
            t.K0s_DeltaE     = (true_K0s->E - best_K0s.E)/true_K0s->E;

            t.Tree->Fill();
        }
    }
}

void TestSigmaPlusFitter::ShowResult()
{

    canvas(GetName())
            << TTree_drawable(t.Tree, "SigmaPlus_DeltaAngle")
            << TTree_drawable(t.Tree, "SigmaPlus_DeltaE")
            << TTree_drawable(t.Tree, "K0s_DeltaAngle")
            << TTree_drawable(t.Tree, "K0s_DeltaE")
            << endc;
}

AUTO_REGISTER_PHYSICS(TestSigmaPlusFitter)