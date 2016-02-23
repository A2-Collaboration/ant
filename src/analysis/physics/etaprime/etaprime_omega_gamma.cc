#include "etaprime_omega_gamma.h"

#include "plot/root_draw.h"
#include "plot/TH1Dcut.h"

#include "utils/particle_tools.h"
#include "utils/matcher.h"
#include "utils/combinatorics.h"

#include "expconfig/ExpConfig.h"

#include "base/Logger.h"
#include "base/std_ext/math.h"

#include <TTree.h>

#include <limits>


using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

const ParticleTypeTree EtapOmegaG::ptreeSignal = ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::EtaPrime_gOmega_ggPi0_4g);
const ParticleTypeTree EtapOmegaG::ptreeReference = ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::EtaPrime_2g);


EtapOmegaG::EtapOmegaG(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    kinfitter_2("kinfitter_2",2),
    kinfitter_4("kinfitter_4",4)
{
    const interval<double> prompt_range{-2.5,1.5};
    promptrandom.AddPromptRange(prompt_range); // slight offset due to CBAvgTime reference
    promptrandom.AddRandomRange({-50,-10});  // just ensure to be way off prompt peak
    promptrandom.AddRandomRange({  10,50});

    promptrandom_tight.AddPromptRange(prompt_range); // slight offset due to CBAvgTime reference
    promptrandom_tight.AddRandomRange({-30,-10});  // just ensure to be way off prompt peak
    promptrandom_tight.AddRandomRange({  10,30});


    h_CommonCuts = HistFac.makeTH1D("Common Cuts", "", "#", BinSettings(15),"h_TotalEvents");
    h_MissedBkg = HistFac.makeTH1D("Missed Background", "", "#", BinSettings(25),"h_MissedBkg");


    const auto setup = ant::ExpConfig::Setup::GetLastFound();
    if(!setup)
        throw runtime_error("EtapOmegaG needs a setup");
    kinfitter_2.LoadSigmaData(setup->GetPhysicsFilesDirectory()+"/FitterSigmas.root");
    kinfitter_4.LoadSigmaData(setup->GetPhysicsFilesDirectory()+"/FitterSigmas.root");

    t.CreateBranches(HistFac.makeTTree("treeCommon"));

    HistogramFactory HistFacNoKinFit("NoKinFit",HistFac);
    HistogramFactory HistFacKinFit("KinFit",HistFac);

    Sig.SetupTrees(HistFacNoKinFit);
    Ref.t.CreateBranches(HistFacNoKinFit.makeTTree("Ref"));

    SigKinFit.SetupTrees(HistFacKinFit);
    RefKinFit.t.CreateBranches(HistFacKinFit.makeTTree("Ref"));

}

void EtapOmegaG::ProcessEvent(const TEvent& event, manager_t&)
{
    // we start with some general candidate handling,
    // later we split into ref/sig analysis according to
    // number of photons

    TEventData& data = *event.Reconstructed;

    h_CommonCuts->Fill("Seen",1.0);

    auto& particletree = event.MCTrue->ParticleTree;

    h_CommonCuts->Fill("MCTrue #eta'", 0); // ensure the bin is there...
    if(particletree) {
        // note: this might also match to g p -> eta' eta' p,
        // but this is kinematically forbidden
        if(utils::ParticleTools::FindParticle(ParticleTypeDatabase::EtaPrime, particletree, 1)) {
            h_CommonCuts->Fill("MCTrue #eta'", 1);
        }
    }

    if(data.Trigger.CBEnergySum<=550)
        return;
    h_CommonCuts->Fill("CBEnergySum>550",1.0);
    t.CBSumE = data.Trigger.CBEnergySum;

    t.CBAvgTime = event.Reconstructed->Trigger.CBTiming;
    if(!isfinite(t.CBAvgTime))
        return;
    h_CommonCuts->Fill("CBAvgTime ok",1.0);

    // use struct to gather particles
    Particles_t particles;

    // identify the proton here as slowest cluster in TAPS
    // using CBAvgTime does not help here, since it's constant
    // over each TAPS clusters
    /// \todo think about using beta here as in EtapProton?
    t.ProtonTime = std_ext::NaN;
    for(const TCandidatePtr& cand : data.Candidates) {
        if(cand->Detector & Detector_t::Type_t::TAPS) {
            if(!isfinite(t.ProtonTime) || t.ProtonTime < cand->Time) {
                t.ProtonTime = cand->Time;
                particles.Proton = make_shared<TParticle>(ParticleTypeDatabase::Proton, cand);
            }
        }
    }

    if(!particles.Proton)
        return;

    h_CommonCuts->Fill("p in TAPS", 1.0);

    // remaining candidates are photons
    const auto nPhotons = data.Candidates.size() - 1;
    if(nPhotons != 2 && nPhotons != 4)
        return;
    h_CommonCuts->Fill("nPhotons==2|4", 1.0);

    particles.PhotonSum.SetPxPyPzE(0,0,0,0);
    t.nPhotonsCB = 0;
    t.nPhotonsTAPS = 0;
    t.CBSumVetoE = 0;
    for(const TCandidatePtr& cand : data.Candidates) {
         if(cand == particles.Proton->Candidate)
             continue;
         if(cand->Detector & Detector_t::Type_t::CB) {
             t.nPhotonsCB++;
             t.CBSumVetoE += cand->VetoEnergy;
         }
         if(cand->Detector & Detector_t::Type_t::TAPS)
             t.nPhotonsTAPS++;
         auto photon = make_shared<TParticle>(ParticleTypeDatabase::Photon, cand);
         particles.PhotonSum += *photon;
         particles.Photons.emplace_back(move(photon));
    }
    assert(particles.Photons.size() == nPhotons);
    assert(nPhotons == t.nPhotonsCB + t.nPhotonsTAPS);

    // don't bother with events where proton coplanarity is not ok
    // we use some rather large window here...
    t.ProtonCopl = std_ext::radian_to_degree(TVector2::Phi_mpi_pi(particles.Proton->Phi() - particles.PhotonSum.Phi() - M_PI ));
    const interval<double> ProtonCopl_cut(-30, 30);
    if(!ProtonCopl_cut.Contains(t.ProtonCopl))
        return;
    h_CommonCuts->Fill("ProtonCopl ok", 1.0);

    // sum up the PID energy
    // (might be different to matched CB/PID Veto energy)
    t.PIDSumE = 0;
    for(const TClusterPtr& cl : data.Clusters) {
        if(cl->DetectorType == Detector_t::Type_t::PID) {
            t.PIDSumE += cl->Energy;
        }
    }



    // select signal or reference according to nPhotons

    t.IsSignal = nPhotons == 4; // else nPhotons==2, so reference
    utils::KinFitter& fitter = t.IsSignal ? kinfitter_4 : kinfitter_2;



    // do some MCTrue identification (if available)
    t.MCTrue = 0; // indicate data by default
    TParticleTree_t ptree_sigref = nullptr; // used by Sig_t::Process to check matching
    if(particletree) {
        auto get_mctrue = [] (ParticleTypeTree expected, TParticleTree_t particletree) -> unsigned {
            if(particletree->IsEqual(expected, utils::ParticleTools::MatchByParticleName))
                return 1;
            unsigned n = 10;
            for(const auto& ptreeBkg : ptreeBackgrounds) {
                if(particletree->IsEqual(ptreeBkg.Tree, utils::ParticleTools::MatchByParticleName))
                    return n;
                n++;
            }
            return 9; // missed background
        };

        // 1=Signal, 2=Reference, 9=MissedBkg, >=10 found in ptreeBackgrounds
        if(t.IsSignal) {
            t.MCTrue = get_mctrue(ptreeSignal, particletree);
            if(t.MCTrue==1)
                ptree_sigref = particletree;
        }
        else {
            t.MCTrue = get_mctrue(ptreeReference, particletree);
            if(t.MCTrue==1) {
                ptree_sigref = particletree;
                t.MCTrue++;
            }
        }

        if(t.MCTrue==9) {
            const auto& decaystr = utils::ParticleTools::GetDecayString(particletree);
            h_MissedBkg->Fill(decaystr.c_str(), 1.0);
        }
    }
    else if(!event.MCTrue->ID.IsInvalid()) {
        // in rare cases, the particletree is not available, although we're running on MCTrue
        // mark this as other MC background
        t.MCTrue = 9;
    }

    Sig.ResetBranches();
    Ref.ResetBranches();
    if(t.IsSignal)
        Sig.Process(particles, ptree_sigref);
    else
        Ref.Process(particles);

    // loop over tagger hits, do KinFit
    bool kinfit_ok = false;
    for(const TTaggerHit& taggerhit : data.TaggerHits) {
        promptrandom.SetTaggerHit(taggerhit.Time - t.CBAvgTime);
        promptrandom_tight.SetTaggerHit(taggerhit.Time - t.CBAvgTime);
        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        // missing mass
        const TLorentzVector beam_target = taggerhit.GetPhotonBeam() + TLorentzVector(0, 0, 0, ParticleTypeDatabase::Proton.Mass());
        t.MissingMass = (beam_target - particles.PhotonSum).M();

        t.TaggW = promptrandom.FillWeight();
        t.TaggW_tight = promptrandom_tight.FillWeight();
        t.TaggE = taggerhit.PhotonEnergy;
        t.TaggT = taggerhit.Time;
        t.TaggCh = taggerhit.Channel;

        // do kinfit
        fitter.SetEgammaBeam(taggerhit.PhotonEnergy);
        fitter.SetProton(particles.Proton);
        fitter.SetPhotons(particles.Photons);
        const auto& fit_result = fitter.DoFit();

        SigKinFit.ResetBranches();
        RefKinFit.ResetBranches();

        t.KinFitChi2 = std_ext::NaN;
        t.KinFitIterations = 0;
        if(fit_result.Status == APLCON::Result_Status_t::Success) {
            kinfit_ok = true;
            t.KinFitChi2 = fit_result.ChiSquare;
            t.KinFitIterations = fit_result.NIterations;
            Particles_t fitted_particles;
            fitted_particles.Proton = fitter.GetFittedProton();
            fitted_particles.Photons = fitter.GetFittedPhotons();
            fitted_particles.PhotonSum.SetPxPyPzE(0,0,0,0);

            for(const TParticlePtr& p : fitted_particles.Photons)
                fitted_particles.PhotonSum += *p;

            if(t.IsSignal)
                SigKinFit.Process(fitted_particles, ptree_sigref);
            else
                RefKinFit.Process(fitted_particles);
        }

        t.Tree->Fill();

        Sig.Fill();
        SigKinFit.Fill();

        Ref.t.Tree->Fill();
        RefKinFit.t.Tree->Fill();

    }

    if(kinfit_ok)
        h_CommonCuts->Fill("KinFit OK", 1.0);
    if(nPhotons==2)
        h_CommonCuts->Fill("nPhotons==2", 1.0);
    if(nPhotons==4)
        h_CommonCuts->Fill("nPhotons==4", 1.0);

}

EtapOmegaG::Sig_t::Sig_t() :
    All(),
    No_Pi0(&ParticleTypeDatabase::Pi0),
    No_Omega(&ParticleTypeDatabase::Omega),
    No_EtaPrime(&ParticleTypeDatabase::EtaPrime)
{
}

void EtapOmegaG::Sig_t::SetupTrees(HistogramFactory HistFac)
{
    All.t.CreateBranches(HistFac.makeTTree("SigAll"));
    No_Pi0.t.CreateBranches(HistFac.makeTTree("SigNoPi0"));
    No_Omega.t.CreateBranches(HistFac.makeTTree("SigNoOmega"));
    No_EtaPrime.t.CreateBranches(HistFac.makeTTree("SigNoEtaPrime"));
}

void EtapOmegaG::Sig_t::Fill()
{
    All.t.Tree->Fill();
    No_Pi0.t.Tree->Fill();
    No_Omega.t.Tree->Fill();
    No_EtaPrime.t.Tree->Fill();
}

void EtapOmegaG::Sig_t::ResetBranches()
{
    All.ResetBranches();
    No_Pi0.ResetBranches();
    No_Omega.ResetBranches();
    No_EtaPrime.ResetBranches();
}

void EtapOmegaG::Sig_t::Process(const Particles_t& particles, TParticleTree_t particletree)
{
    All.Process(particles, particletree);
    No_Pi0.Process(particles, particletree);
    No_Omega.Process(particles, particletree);
    No_EtaPrime.Process(particles, particletree);
}

EtapOmegaG::Sig_t::Fit_t::Fit_t(const ParticleTypeDatabase::Type* typeptr) :
    treefitter{MakeFitter(typeptr)},
    fitted_EtaPrime(treefitter.GetTreeNode(ParticleTypeDatabase::EtaPrime)),
    fitted_Omega(treefitter.GetTreeNode(ParticleTypeDatabase::Omega)),
    fitted_Pi0(treefitter.GetTreeNode(ParticleTypeDatabase::Pi0))
{
    t.ggg().resize(4);
    t.gg_gg1().resize(3);
    t.gg_gg2().resize(3);

    const auto setup = ant::ExpConfig::Setup::GetLastFound();
    if(!setup)
        throw runtime_error("EtapOmegaG needs a setup");
    treefitter.LoadSigmaData(setup->GetPhysicsFilesDirectory()+"/FitterSigmas.root");

    // search dependent gammas and remember the tree nodes in the fitter

    auto find_photons = [] (utils::TreeFitter::tree_t fitted) {
        std::vector<utils::TreeFitter::tree_t> photons;
        for(const auto& d : fitted->Daughters())
            if(d->Get().TypeTree->Get() == ParticleTypeDatabase::Photon)
                photons.emplace_back(d);
        return photons;
    };

    fitted_g_EtaPrime = find_photons(fitted_EtaPrime).at(0);
    fitted_g_Omega = find_photons(fitted_Omega).at(0);
    fitted_g1_Pi0 = find_photons(fitted_Pi0).at(0);
    fitted_g2_Pi0 = find_photons(fitted_Pi0).at(1);
}

utils::TreeFitter EtapOmegaG::Sig_t::Fit_t::MakeFitter(const ParticleTypeDatabase::Type* typeptr)
{
    if(typeptr == nullptr)
        return {"sig_treefitter_All", utils::ParticleTools::GetProducedParticle(EtapOmegaG::ptreeSignal)};
    else
        return {"sig_treefitter_No"+typeptr->Name(),
                    utils::ParticleTools::GetProducedParticle(EtapOmegaG::ptreeSignal),
                    [&typeptr] (ParticleTypeTree tree) { return tree->Get() == *typeptr; }
        };
}

void EtapOmegaG::Sig_t::Fit_t::ResetBranches()
{
    t.TreeFitChi2 = std_ext::NaN;
    t.TreeFitIterations = 0;

    std::fill(t.ggg().begin(), t.ggg().end(), std_ext::NaN);
    std::fill(t.gg_gg1().begin(), t.gg_gg1().end(), std_ext::NaN);
    std::fill(t.gg_gg2().begin(), t.gg_gg2().end(), std_ext::NaN);

    t.IM_EtaPrime_fitted  = std_ext::NaN;
    t.IM_Omega_fitted = std_ext::NaN;
    t.IM_Pi0_fitted = std_ext::NaN;

    t.IM_EtaPrime_best = std_ext::NaN;
    t.IM_Omega_best = std_ext::NaN;
    t.IM_Pi0_best = std_ext::NaN;

    t.Bachelor_best_best = std_ext::NaN;
    t.Bachelor_best_fit = std_ext::NaN;
    t.Bachelor_fit_best = std_ext::NaN;
    t.Bachelor_fit_fit = std_ext::NaN;

    t.MCTrueMatch = 0;
}

void EtapOmegaG::Sig_t::Fit_t::Process(const EtapOmegaG::Particles_t& particles, TParticleTree_t ptree_sigref)
{
    ResetBranches();

    assert(particles.Photons.size() == 4);

    //  ggg combinatorics
    auto it_ggg = t.ggg().begin();
    for( auto comb = utils::makeCombination(particles.Photons,3); !comb.Done(); ++comb ) {
        *it_ggg = (*comb.at(0) + *comb.at(1) + *comb.at(2)).M();
        ++it_ggg;
    }

    // gg/gg "Goldhaber" combinatorics
    const auto goldhaber_comb = vector<vector<unsigned>>({{0,1,2,3},{0,2,1,3},{0,3,1,2}});
    auto it_gg_gg1 = t.gg_gg1().begin();
    auto it_gg_gg2 = t.gg_gg2().begin();
    for(auto i : goldhaber_comb) {
        const auto& p = particles.Photons;
        *it_gg_gg1 = (*p[i[0]] + *p[i[1]]).M();
        *it_gg_gg2 = (*p[i[2]] + *p[i[3]]).M();
        ++it_gg_gg1;
        ++it_gg_gg2;
    }


    // g_Omega to check against MCTrue
    TParticlePtr g_Omega_best;
    // the EtaPrime bachelor photon is most important to us...
    TParticlePtr g_EtaPrime_best;
    TLorentzVector g_EtaPrime_fitted;
    const TLorentzVector& EtaPrime_best = particles.PhotonSum; // does not change with permutation
    t.IM_EtaPrime_best = EtaPrime_best.M();
    TLorentzVector EtaPrime_fitted;

    treefitter.SetLeaves(particles.Photons);
    APLCON::Result_t r;

    while(treefitter.NextFit(r)) {
        if(r.Status != APLCON::Result_Status_t::Success)
            continue;
        if(isfinite(t.TreeFitChi2) && r.ChiSquare>t.TreeFitChi2)
            continue;
        // found fit with better chi2
        t.TreeFitChi2 = r.ChiSquare;
        t.TreeFitIterations = r.NIterations;

        EtaPrime_fitted = fitted_EtaPrime->Get().Particle;

        // IM fitted expected to be delta peaks since they were fitted...
        t.IM_EtaPrime_fitted = EtaPrime_fitted.M();
        t.IM_Omega_fitted = fitted_Omega->Get().Particle.M();
        t.IM_Pi0_fitted = fitted_Pi0->Get().Particle.M();

        // have a look at the assigned gammas to Pi0/Omega
        const TLorentzVector& Pi0_best = *fitted_g1_Pi0->Get().SetParticle + *fitted_g2_Pi0->Get().SetParticle;
        t.IM_Pi0_best = Pi0_best.M();

        g_Omega_best = fitted_g_Omega->Get().SetParticle;
        const TLorentzVector& Omega_best = *g_Omega_best + Pi0_best;
        t.IM_Omega_best = Omega_best.M();

        // have a look at the EtaPrime bachelor photon
        g_EtaPrime_fitted = fitted_g_EtaPrime->Get().Particle;
        g_EtaPrime_best = fitted_g_EtaPrime->Get().SetParticle;
    }

    if(isfinite(t.TreeFitChi2)) {

        // do that Bachelor photon boosting
        auto do_boost = [] (const TLorentzVector& bachelor, const TLorentzVector& etaprime) {
            TLorentzVector boosted(bachelor);
            boosted.Boost(-etaprime.BoostVector());
            return boosted;
        };

        t.Bachelor_best_best = do_boost(*g_EtaPrime_best, EtaPrime_best).E();
        t.Bachelor_best_fit  = do_boost(*g_EtaPrime_best, EtaPrime_fitted).E();
        t.Bachelor_fit_best  = do_boost(g_EtaPrime_fitted, EtaPrime_best).E();
        t.Bachelor_fit_fit   = do_boost(g_EtaPrime_fitted, EtaPrime_fitted).E();


        // check MC matching
        if(ptree_sigref) {

            auto true_photons = utils::ParticleTools::FindParticles(ParticleTypeDatabase::Photon, ptree_sigref);
            assert(true_photons.size() == 4);
            auto match_bycandidate = [] (const TParticlePtr& mctrue, const TParticlePtr& recon) {
                return mctrue->Angle(*recon->Candidate); // TCandidate converts into TVector3
            };
            auto matched = utils::match1to1(true_photons, particles.Photons,
                                            match_bycandidate,IntervalD(0.0, 15.0*TMath::DegToRad()));
            if(matched.size() == 4) {
                // do that tedious photon determination (rewriting the matcher somehow would be nice....)
                auto select_daughter = [] (TParticleTree_t tree, const ParticleTypeDatabase::Type& type) {
                    auto d = tree->Daughters().front()->Get()->Type() == type ?
                                 tree->Daughters().front() : tree->Daughters().back();
                    assert(d->Get()->Type() == type);
                    return d;
                };

                auto etap = select_daughter(ptree_sigref, ParticleTypeDatabase::EtaPrime);
                auto g_Etap = select_daughter(etap, ParticleTypeDatabase::Photon);
                auto omega = select_daughter(etap, ParticleTypeDatabase::Omega);
                auto g_Omega = select_daughter(omega, ParticleTypeDatabase::Photon);

                auto g_Etap_matched = utils::FindMatched(matched, g_Etap->Get());
                auto g_Omega_matched = utils::FindMatched(matched, g_Omega->Get());
                if(g_Etap_matched == g_EtaPrime_best)
                    t.MCTrueMatch += 1;
                if(g_Omega_matched == g_Omega_best)
                    t.MCTrueMatch += 2;

            }
        }
    }


}

void EtapOmegaG::Ref_t::ResetBranches()
{
    t.IM_2g = std_ext::NaN;
}

void EtapOmegaG::Ref_t::Process(const EtapOmegaG::Particles_t& particles)
{
    assert(particles.Photons.size() == 2);
    t.IM_2g = particles.PhotonSum.M();
}

void EtapOmegaG::ShowResult()
{
    canvas("Overview") << h_CommonCuts << h_MissedBkg << endc;
}


const std::vector<EtapOmegaG::Background_t> EtapOmegaG::ptreeBackgrounds = {
    {"1Pi0", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0_2g)},
    {"2Pi0", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::TwoPi0_4g)},
    {"Pi0Eta", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0Eta_4g)},
    {"3Pi0", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::ThreePi0_6g)},
    {"OmegaPi0g", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Omega_gPi0_3g)},
    {"OmegaPi0PiPPiM", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Omega_Pi0PiPPiM_2g)},
    {"EtaP2Pi0Eta", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::EtaPrime_2Pi0Eta_6g)},
    {"2Pi0Dalitz", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::TwoPi0_2ggEpEm)},
    {"3Pi0Dalitz", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::ThreePi0_4ggEpEm)},
    {"1Eta", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Eta_2g)},
};

AUTO_REGISTER_PHYSICS(EtapOmegaG)
