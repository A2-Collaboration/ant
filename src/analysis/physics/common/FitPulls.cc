#include "FitPulls.h"
#include "base/std_ext/vector.h"
#include "TH1D.h"

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace std;

const vector<ParticleTypeTree> FitPulls::channels = {
    ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0_2g),
    ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::TwoPi0_4g),
    ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::ThreePi0_6g),
    ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Omega_gPi0_3g),
};

FitPulls::FitPulls(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    opt_save_after_cut(opts->Get<bool>("SaveAfterCut", false)),
    opt_save_only(opts->Get<bool>("SaveOnly", false))
{
    promptrandom.AddPromptRange({-2.5,1.5}); // slight offset due to CBAvgTime reference
    promptrandom.AddRandomRange({-50,-10});  // just ensure to be way off prompt peak
    promptrandom.AddRandomRange({  10,50});

    h_protoncopl = HistFac.makeTH1D("p Coplanarity","#delta#phi / degree","",BinSettings(400,-180,180),"h_protoncopl");
    h_taggtime = HistFac.makeTH1D("Tagged Time","t / ns", "", BinSettings(300,-60,60), "h_taggtime");


    // create fitter for each channel
    auto count_leaves = [] (const ParticleTypeTree& tree) {
        unsigned leaves = 0;
        tree->Map_nodes(
                    [&leaves] (const ParticleTypeTree& n) { if(n->IsLeaf()) leaves++; }
        );
        return leaves;
    };

    auto uncertainty_model = make_shared<utils::UncertaintyModels::Optimized_Oli1>();

    unsigned n_ch = 0;
    for(auto& channel : channels) {
        const auto nParticles = count_leaves(channel);
        treefitters[nParticles].emplace_back(
                    std_ext::formatter() << "treefitter_" << n_ch,
                    channel,
                    nParticles-1, // assume one particle to be proton
                    uncertainty_model
                    );
        n_ch++;
    }
}

void FitPulls::findProton(const TCandidateList& cands,
                          const TTaggerHit& taggerhit,
                          TParticlePtr& proton,
                          TParticleList& photons,
                          LorentzVec& photon_sum)
{
    // find the candidate which gives the closest missing mass to proton
    auto min_mm = std::numeric_limits<double>::infinity();
    for(auto& cand_proton : cands.get_iter()) {
        LorentzVec photon_sum_tmp(0,0,0,0);
        for(auto& cand : cands.get_iter()) {
            if(cand == cand_proton)
                continue;
            photon_sum_tmp += TParticle(ParticleTypeDatabase::Photon, cand.get_ptr());
        }
        const LorentzVec beam_target = taggerhit.GetPhotonBeam() + LorentzVec(0, 0, 0, ParticleTypeDatabase::Proton.Mass());
        const LorentzVec missing = beam_target - photon_sum_tmp;
        const double mm = abs(missing.M() - ParticleTypeDatabase::Proton.Mass());
        if(mm < min_mm) {
            min_mm = mm;
            proton = make_shared<TParticle>(ParticleTypeDatabase::Proton, cand_proton.get_ptr());
            photon_sum = photon_sum_tmp;
        }
    }

    // fill the photons
    for(auto& cand_photon : cands.get_iter()) {
        auto cand_ptr = cand_photon.get_ptr();
        if(cand_ptr == proton->Candidate)
            continue;
        photons.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::Photon, cand_ptr));
    }
}

void FitPulls::ProcessEvent(const TEvent& event, manager_t& manager)
{
    const auto& cands = event.Reconstructed().Candidates;

    auto fitters = treefitters.find(cands.size());
    if(fitters == treefitters.end())
        return;

    // check that all candidates are clean
    for(const TCandidate& cand : cands) {
        if(auto calocluster = cand.FindCaloCluster()) {
            if(calocluster->HasFlag(TCluster::Flags_t::TouchesHole))
                return;
        }
        else
            return;
    }

    if(opt_save_after_cut) {
        manager.SaveEvent();
        if(opt_save_only)
            return;
    }

    for(const TTaggerHit& taggerhit : event.Reconstructed().TaggerHits) {
        const auto taggtime = taggerhit.Time - event.Reconstructed().Trigger.CBTiming;
        promptrandom.SetTaggerHit(taggtime);
        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        h_taggtime->Fill(taggtime);

        // select the proton as the particle with best missing mass
        TParticlePtr  proton;
        TParticleList photons;
        LorentzVec photon_sum;
        findProton(cands, taggerhit, proton, photons, photon_sum);

        // check coplanarity
        const double d_phi = std_ext::radian_to_degree(vec2::Phi_mpi_pi(proton->Phi()-photon_sum.Phi() - M_PI ));
        const interval<double> copl_cut{-20, 20};
        if(!copl_cut.Contains(d_phi))
            continue;

        h_protoncopl->Fill(d_phi, promptrandom.FillWeight());

        // do all the fits
        for(utils::TreeFitter& fitter : fitters->second) {
            fitter.SetEgammaBeam(taggerhit.PhotonEnergy);
            fitter.SetPhotons(photons);
            fitter.SetProton(proton);

            APLCON::Result_t result;
            while(fitter.NextFit(result)) {

            }
        }
    }
}

void FitPulls::ShowResult()
{

}



AUTO_REGISTER_PHYSICS(FitPulls)