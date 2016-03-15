#include "etaprime_omega_gamma.h"

#include "plot/root_draw.h"

#include "utils/particle_tools.h"
#include "utils/matcher.h"
#include "utils/combinatorics.h"

#include "expconfig/ExpConfig.h"

#include "base/Logger.h"
#include "base/std_ext/math.h"
#include "base/std_ext/vector.h"

#include <TTree.h>

#include <limits>


using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

const ParticleTypeTree EtapOmegaG::ptreeSignal = ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::EtaPrime_gOmega_ggPi0_4g);
const ParticleTypeTree EtapOmegaG::ptreeReference = ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::EtaPrime_2g);


APLCON::Fit_Settings_t EtapOmegaG::MakeFitSettings(unsigned max_iterations)
{
    auto settings = APLCON::Fit_Settings_t::Default;
    settings.MaxIterations = max_iterations;
    return settings;
}

EtapOmegaG::EtapOmegaG(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    kinfitter_2("kinfitter_2", 2, MakeFitSettings(25)),
    kinfitter_4("kinfitter_4", 4, MakeFitSettings(25))
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

    t.CreateBranches(HistFac.makeTTree("treeCommon"));

    Sig.SetupTrees(HistFac);
    Ref.t.CreateBranches(HistFac.makeTTree("Ref"));

    const auto setup = ant::ExpConfig::Setup::GetLastFound();
    if(!setup)
        throw runtime_error("EtapOmegaG needs a setup");
    kinfitter_2.LoadSigmaData(setup->GetPhysicsFilesDirectory()+"/FitterSigmas.root");
    kinfitter_4.LoadSigmaData(setup->GetPhysicsFilesDirectory()+"/FitterSigmas.root");

}

void EtapOmegaG::ProcessEvent(const TEvent& event, manager_t&)
{
    // we start with some general candidate handling,
    // later we split into ref/sig analysis according to
    // number of photons

    const TEventData& data = event.Reconstructed();

    h_CommonCuts->Fill("Seen",1.0);

    auto& particletree = event.MCTrue().ParticleTree;

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

    t.CBAvgTime = event.Reconstructed().Trigger.CBTiming;
    if(!isfinite(t.CBAvgTime))
        return;
    h_CommonCuts->Fill("CBAvgTime ok",1.0);

    // use struct to gather particles

    // identify the proton here as slowest cluster in TAPS
    // using CBAvgTime does not help here, since it's constant
    // over each TAPS clusters
    /// \todo think about using beta here as in EtapProton?
    t.ProtonTime = std_ext::NaN;
    TParticlePtr proton;
    for(const TCandidatePtr& cand : data.Candidates) {
        if(cand->Detector & Detector_t::Type_t::TAPS) {
            if(!isfinite(t.ProtonTime) || t.ProtonTime < cand->Time) {
                t.ProtonTime = cand->Time;
                proton = make_shared<TParticle>(ParticleTypeDatabase::Proton, cand);
            }
        }
    }

    if(!proton)
        return;
    h_CommonCuts->Fill("p in TAPS", 1.0);

    // remaining candidates are photons
    const auto nPhotons = data.Candidates.size() - 1;
    if(nPhotons != 2 && nPhotons != 4)
        return;
    h_CommonCuts->Fill("nPhotons==2|4", 1.0);

    Particles_t particles;
    TLorentzVector photon_sum(0,0,0,0);
    t.nPhotonsCB = 0;
    t.nPhotonsTAPS = 0;
    t.CBSumVetoE = 0;
    for(const TCandidatePtr& cand : data.Candidates) {
         if(cand == proton->Candidate)
             continue;
         if(cand->Detector & Detector_t::Type_t::CB) {
             t.nPhotonsCB++;
             t.CBSumVetoE += cand->VetoEnergy;
         }
         if(cand->Detector & Detector_t::Type_t::TAPS)
             t.nPhotonsTAPS++;
         auto photon = make_shared<TParticle>(ParticleTypeDatabase::Photon, cand);
         photon_sum += *photon;
         particles.Photons.emplace_back(move(photon));
    }
    assert(particles.Photons.size() == nPhotons);
    assert(nPhotons == t.nPhotonsCB + t.nPhotonsTAPS);

    // don't bother with events where proton coplanarity is not ok
    // we use some rather large window here...
    t.ProtonCopl = std_ext::radian_to_degree(TVector2::Phi_mpi_pi(proton->Phi() - photon_sum.Phi() - M_PI ));
    const interval<double> ProtonCopl_cut(-35, 35);
    if(!ProtonCopl_cut.Contains(t.ProtonCopl))
        return;
    h_CommonCuts->Fill("ProtonCopl ok", 1.0);

    // sum up the PID energy
    // (might be different to matched CB/PID Veto energy)
    t.PIDSumE = 0;
    for(const TCluster& cl : data.Clusters) {
        if(cl.DetectorType == Detector_t::Type_t::PID) {
            t.PIDSumE += cl.Energy;
        }
    }



    // select signal or reference according to nPhotons
    t.IsSignal = nPhotons == 4; // else nPhotons==2, so reference

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

        // fill histogram of missed bkgs
        if(t.MCTrue==9) {
            const auto& decaystr = utils::ParticleTools::GetDecayString(particletree);
            h_MissedBkg->Fill(decaystr.c_str(), 1.0);
        }
    }
    else if(!event.MCTrue().ID.IsInvalid()) {
        // in rare cases, the particletree is not available, although we're running on MCTrue
        // mark this as other MC background
        t.MCTrue = 9;
    }

    auto& kinfitter = nPhotons==4 ? kinfitter_4 : kinfitter_2;

    // loop over tagger hits, delegate to Ref/Sig
    for(const TTaggerHit& taggerhit : data.TaggerHits) {
        promptrandom.SetTaggerHit(taggerhit.Time - t.CBAvgTime);
        promptrandom_tight.SetTaggerHit(taggerhit.Time - t.CBAvgTime);
        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        // missing mass
        const TLorentzVector beam_target = taggerhit.GetPhotonBeam() + TLorentzVector(0, 0, 0, ParticleTypeDatabase::Proton.Mass());
        t.MissingMass = (beam_target - photon_sum).M();

        t.TaggW = promptrandom.FillWeight();
        t.TaggW_tight = promptrandom_tight.FillWeight();
        t.TaggE = taggerhit.PhotonEnergy;
        t.TaggT = taggerhit.Time;
        t.TaggCh = taggerhit.Channel;

        Sig.ResetBranches();
        Ref.ResetBranches();

        // do kinfit
        kinfitter.SetEgammaBeam(taggerhit.PhotonEnergy);
        kinfitter.SetProton(proton);
        kinfitter.SetPhotons(particles.Photons);
        const auto& fit_result = kinfitter.DoFit();

        t.KinFitChi2 = std_ext::NaN;
        t.KinFitProb = std_ext::NaN;
        t.KinFitIterations = 0;

        if(fit_result.Status == APLCON::Result_Status_t::Success) {

            t.KinFitChi2 = fit_result.ChiSquare;
            t.KinFitProb = fit_result.Probability;
            t.KinFitIterations = fit_result.NIterations;

            Particles_t fitted_particles;
            fitted_particles.Photons = kinfitter.GetFittedPhotons();
            fitted_particles.PhotonSum.SetPxPyPzE(0,0,0,0);
            for(const auto& photon : fitted_particles.Photons)
                fitted_particles.PhotonSum += *photon;

            if(t.IsSignal)
                Sig.Process(fitted_particles, ptree_sigref);
            else
                Ref.Process(fitted_particles);
        }

        t.Tree->Fill();
        Sig.Fill();
        Ref.t.Tree->Fill();
    }

    if(nPhotons==2)
        h_CommonCuts->Fill("nPhotons==2", 1.0);
    if(nPhotons==4)
        h_CommonCuts->Fill("nPhotons==4", 1.0);

}

EtapOmegaG::Sig_t::Sig_t() :
    treefitter_Pi0Pi0("treefit_Pi0Pi0",
                      ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::TwoPi0_4g), {},
                      MakeFitSettings(20)
                      ),
    treefitter_Pi0Eta("treefit_Pi0Eta",
                      ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0Eta_4g), {},
                      MakeFitSettings(15)
                      )
{
    const auto setup = ant::ExpConfig::Setup::GetLastFound();
    if(!setup)
        throw runtime_error("EtapOmegaG needs a setup");
    treefitter_Pi0Pi0.LoadSigmaData(setup->GetPhysicsFilesDirectory()+"/FitterSigmas.root");
    treefitter_Pi0Eta.LoadSigmaData(setup->GetPhysicsFilesDirectory()+"/FitterSigmas.root");
}

void EtapOmegaG::Sig_t::SetupTrees(HistogramFactory HistFac)
{
    t.CreateBranches(HistFac.makeTTree("SigShared"));
    OmegaPi0.t.CreateBranches(HistFac.makeTTree("SigOmegaPi0"));
    Pi0.t.CreateBranches(HistFac.makeTTree("SigPi0"));
}

void EtapOmegaG::Sig_t::Fill()
{
    t.Tree->Fill();
    OmegaPi0.t.Tree->Fill();
    Pi0.t.Tree->Fill();
}

void EtapOmegaG::Sig_t::ResetBranches()
{
    t.Reset();
    OmegaPi0.t.Reset();
    Pi0.t.Reset();
}

void EtapOmegaG::Sig_t::SharedTree_t::Reset()
{
    std::fill(ggg().begin(), ggg().end(), std_ext::NaN);
    std::fill(gg_gg1().begin(), gg_gg1().end(), std_ext::NaN);
    std::fill(gg_gg2().begin(), gg_gg2().end(), std_ext::NaN);

    AntiPi0FitChi2 = std_ext::NaN;
    AntiPi0FitProb = std_ext::NaN;
    AntiPi0FitIterations = 0;

    AntiEtaFitChi2 = std_ext::NaN;
    AntiEtaFitProb = std_ext::NaN;
    AntiEtaFitIterations = 0;
}

void EtapOmegaG::Sig_t::Process(const Particles_t& particles, TParticleTree_t ptree_sigref)
{
    DoPhotonCombinatorics(particles.Photons);
    DoAntiPi0Eta(particles.Photons);
    OmegaPi0.Process(particles, ptree_sigref);
    Pi0.Process(particles, ptree_sigref);
}

void EtapOmegaG::Sig_t::DoPhotonCombinatorics(TParticleList photons)
{

    //  ggg combinatorics
    auto it_ggg = t.ggg().begin();
    for( auto comb = utils::makeCombination(photons,3); !comb.Done(); ++comb ) {
        *it_ggg = (*comb.at(0) + *comb.at(1) + *comb.at(2)).M();
        ++it_ggg;
    }

    // gg/gg "Goldhaber" combinatorics
    const auto goldhaber_comb = vector<vector<unsigned>>({{0,1,2,3},{0,2,1,3},{0,3,1,2}});
    auto it_gg_gg1 = t.gg_gg1().begin();
    auto it_gg_gg2 = t.gg_gg2().begin();
    for(auto i : goldhaber_comb) {
        const auto& p = photons;
        *it_gg_gg1 = (*p[i[0]] + *p[i[1]]).M();
        *it_gg_gg2 = (*p[i[2]] + *p[i[3]]).M();
        ++it_gg_gg1;
        ++it_gg_gg2;
    }
}

void EtapOmegaG::Sig_t::DoAntiPi0Eta(TParticleList photons)
{
    APLCON::Result_t r;

    treefitter_Pi0Pi0.SetLeaves(photons);
    while(treefitter_Pi0Pi0.NextFit(r)) {
        if(r.Status != APLCON::Result_Status_t::Success)
            continue;
        if(isfinite(t.AntiPi0FitChi2) && r.ChiSquare>t.AntiPi0FitChi2)
            continue;
        // found fit with better chi2
        t.AntiPi0FitChi2 = r.ChiSquare;
        t.AntiPi0FitProb = r.Probability;
        t.AntiPi0FitIterations = r.NIterations;
    }

    treefitter_Pi0Eta.SetLeaves(photons);
    while(treefitter_Pi0Eta.NextFit(r)) {
        if(r.Status != APLCON::Result_Status_t::Success)
            continue;
        if(isfinite(t.AntiEtaFitChi2) && r.ChiSquare>t.AntiEtaFitChi2)
            continue;
        // found fit with better chi2
        t.AntiEtaFitChi2 = r.ChiSquare;
        t.AntiEtaFitProb = r.Probability;
        t.AntiEtaFitIterations = r.NIterations;
    }
}

EtapOmegaG::Sig_t::Fit_t::Fit_t(utils::TreeFitter fitter) :
    treefitter(move(fitter)),
    fitted_Pi0(treefitter.GetTreeNode(ParticleTypeDatabase::Pi0)),
    fitted_Omega(treefitter.GetTreeNode(ParticleTypeDatabase::Omega))
{
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

    fitted_g1_Pi0 = find_photons(fitted_Pi0).at(0);
    fitted_g2_Pi0 = find_photons(fitted_Pi0).at(1);

    if(fitted_Omega)
        fitted_g_Omega = find_photons(fitted_Omega).at(0);
}

utils::TreeFitter EtapOmegaG::Sig_t::Fit_t::Make(const ParticleTypeDatabase::Type& subtree)
{
    auto EtaPrime = utils::ParticleTools::GetProducedParticle(EtapOmegaG::ptreeSignal);

    auto select_daughter = [] (ParticleTypeTree tree, const ParticleTypeDatabase::Type& type) {
        auto d = tree->Daughters().front()->Get() == type ?
                     tree->Daughters().front() : tree->Daughters().back();
        assert(d->Get() == type);
        return d;
    };

    auto Omega = select_daughter(EtaPrime, ParticleTypeDatabase::Omega);
    auto Pi0 = select_daughter(Omega, ParticleTypeDatabase::Pi0);

    return {
        "sig_treefitter_"+subtree.Name(),
        subtree==Omega->Get() ? Omega : Pi0,
        {},
        MakeFitSettings(15)
    };
}




void EtapOmegaG::Sig_t::Fit_t::BaseTree_t::Reset()
{
    TreeFitChi2 = std_ext::NaN;
    TreeFitProb = std_ext::NaN;
    TreeFitIterations = 0;

    ClusterShape_g1_Pi0 = std_ext::NaN;
    ClusterShape_g2_Pi0 = std_ext::NaN;
    ClusterShape_g_Omega = std_ext::NaN;
    ClusterShape_g_EtaPrime = std_ext::NaN;

    IM_Pi0_best = std_ext::NaN;
    IM_Pi0_fitted = std_ext::NaN;
    IM_Pi0gg = std_ext::NaN;
    IM_gg = std_ext::NaN;

    MCTrueMatch = 0;
}

EtapOmegaG::Sig_t::Pi0_t::Pi0_t() :
    Fit_t(Fit_t::Make(ParticleTypeDatabase::Pi0))
{

}

void EtapOmegaG::Sig_t::Pi0_t::BaseTree_t::Reset()
{
    Fit_t::BaseTree_t::Reset();
    std::fill(IM_Pi0g().begin(), IM_Pi0g().end(), std_ext::NaN);
}

void EtapOmegaG::Sig_t::Pi0_t::Process(const EtapOmegaG::Particles_t& particles,
                                       TParticleTree_t ptree_sigref)
{
    assert(particles.Photons.size() == 4);

    // sum of photons should give the EtaPrime
    t.IM_Pi0gg = particles.PhotonSum.M();

    // for MCtrue identification
    TParticlePtr g1_Pi0_best;
    TParticlePtr g2_Pi0_best;

    for(auto comb=utils::makeCombination(particles.Photons, 2); !comb.Done(); ++comb) {

        // do treefit (build photons from iterator ctor)
        treefitter.SetLeaves({comb.begin(), comb.end()});

        APLCON::Result_t r;

        while(treefitter.NextFit(r)) {
            if(r.Status != APLCON::Result_Status_t::Success)
                continue;
            if(isfinite(t.TreeFitChi2) && r.ChiSquare>t.TreeFitChi2)
                continue;
            // found fit with better chi2
            t.TreeFitChi2 = r.ChiSquare;
            t.TreeFitProb = r.Probability;
            t.TreeFitIterations = r.NIterations;

            // IM fitted expected to be delta peaks since they were fitted...
            const TLorentzVector& pi0_fitted = fitted_Pi0->Get().LVSum;
            t.IM_Pi0_fitted = pi0_fitted.M();

            // have a look at the assigned gammas to Pi0/Omega
            g1_Pi0_best = fitted_g1_Pi0->Get().Leave->Particle;
            g2_Pi0_best = fitted_g2_Pi0->Get().Leave->Particle;
            t.ClusterShape_g1_Pi0 = clustertools.LateralMoment(*g1_Pi0_best->Candidate->FindCaloCluster());
            t.ClusterShape_g2_Pi0 = clustertools.LateralMoment(*g2_Pi0_best->Candidate->FindCaloCluster());

            const TLorentzVector& Pi0_best = *g1_Pi0_best + *g2_Pi0_best;
            t.IM_Pi0_best = Pi0_best.M();

            // there are two photon combinations possible
            // for the omega
            assert(t.IM_Pi0g().size() == 2);
            const TParticlePtr& g1 = *comb.begin_not();
            const TParticlePtr& g2 = *std::next(comb.begin_not());
            t.IM_gg = (*g1 + *g2).M();
            t.IM_Pi0g().front() = (Pi0_best + *g1).M();
            t.IM_Pi0g().back()  = (Pi0_best + *g2).M();
            t.ClusterShape_g_EtaPrime = clustertools.LateralMoment(*g1->Candidate->FindCaloCluster());
            t.ClusterShape_g_Omega    = clustertools.LateralMoment(*g2->Candidate->FindCaloCluster());
            if(t.IM_Pi0g().front() > t.IM_Pi0g().back()) {
                // first IM is higher, then swap also cluster shapes
                // we assume that the photon with the lower IM corresponds to
                // the EtaPrime bachelor photon
                std::swap(t.IM_Pi0g().front(), t.IM_Pi0g().back());
                std::swap(t.ClusterShape_g_EtaPrime(), t.ClusterShape_g_Omega());
            }
        }
    }

    // there was at least one successful fit
    if(isfinite(t.TreeFitChi2)) {

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
                // find the two photons of the pi0
                TParticleList pi0_photons;
                ptree_sigref->Map_nodes([&pi0_photons] (const TParticleTree_t& t) {
                    const auto& parent = t->GetParent();
                    if(!parent)
                        return;
                    if(parent->Get()->Type() == ParticleTypeDatabase::Pi0) {
                        pi0_photons.push_back(t->Get());
                    }
                });
                TParticleList g_pi0_matched{
                    utils::FindMatched(matched, pi0_photons.front()),
                    utils::FindMatched(matched, pi0_photons.back())
                };

                if(std_ext::contains(g_pi0_matched, g1_Pi0_best))
                    t.MCTrueMatch += 1;
                if(std_ext::contains(g_pi0_matched, g2_Pi0_best))
                    t.MCTrueMatch += 2;
            }
        }

    }

}


EtapOmegaG::Sig_t::OmegaPi0_t::OmegaPi0_t() :
    Fit_t(Fit_t::Make(ParticleTypeDatabase::Omega))
{

}

void EtapOmegaG::Sig_t::OmegaPi0_t::BaseTree_t::Reset()
{
    Fit_t::BaseTree_t::Reset();
    IM_Pi0g_fitted = std_ext::NaN;
    IM_Pi0g_best = std_ext::NaN;
    Bachelor_E = std_ext::NaN;
}


void EtapOmegaG::Sig_t::OmegaPi0_t::Process(const EtapOmegaG::Particles_t& particles, TParticleTree_t ptree_sigref)
{

    assert(particles.Photons.size() == 4);

    // sum of photons should give the EtaPrime
    const TLorentzVector& EtaPrime = particles.PhotonSum;
    t.IM_Pi0gg = EtaPrime.M();

    // g_Omega to check against MCTrue
    TParticlePtr g_Omega_best;
    // the EtaPrime bachelor photon is most important to us...
    TParticlePtr g_EtaPrime_best;

    for(auto comb=utils::makeCombination(particles.Photons, 3); !comb.Done(); ++comb) {

        // do treefit (build photons from iterator ctor)
        treefitter.SetLeaves({comb.begin(), comb.end()});

        APLCON::Result_t r;

        while(treefitter.NextFit(r)) {
            if(r.Status != APLCON::Result_Status_t::Success)
                continue;
            if(isfinite(t.TreeFitChi2) && r.ChiSquare>t.TreeFitChi2)
                continue;
            // found fit with better chi2
            t.TreeFitChi2 = r.ChiSquare;
            t.TreeFitProb = r.Probability;
            t.TreeFitIterations = r.NIterations;

            // IM fitted expected to be delta peaks since they were fitted...
            t.IM_Pi0g_fitted = fitted_Omega->Get().LVSum.M();
            t.IM_Pi0_fitted = fitted_Pi0->Get().LVSum.M();

            // have a look at the assigned gammas to Pi0/Omega
            const auto& g1_Pi0_best = fitted_g1_Pi0->Get().Leave->Particle;
            const auto& g2_Pi0_best = fitted_g2_Pi0->Get().Leave->Particle;
            t.ClusterShape_g1_Pi0 = clustertools.LateralMoment(*g1_Pi0_best->Candidate->FindCaloCluster());
            t.ClusterShape_g2_Pi0 = clustertools.LateralMoment(*g2_Pi0_best->Candidate->FindCaloCluster());

            const TLorentzVector& Pi0_best = *g1_Pi0_best + *g2_Pi0_best;
            t.IM_Pi0_best = Pi0_best.M();

            g_Omega_best = fitted_g_Omega->Get().Leave->Particle;
            const TLorentzVector& Omega_best = *g_Omega_best + Pi0_best;
            t.IM_Pi0g_best = Omega_best.M();

            // have a look at the EtaPrime bachelor photon
            // the element NOT in the combination is the Bachelor photon
            g_EtaPrime_best = *comb.begin_not();
            t.IM_gg = (*g_EtaPrime_best + *g_Omega_best).M();

            t.ClusterShape_g_Omega = clustertools.LateralMoment(*g_Omega_best->Candidate->FindCaloCluster());
            t.ClusterShape_g_EtaPrime = clustertools.LateralMoment(*g_EtaPrime_best->Candidate->FindCaloCluster());
        }
    }


    // there was at least one successful fit
    if(isfinite(t.TreeFitChi2)) {

        // do that Bachelor photon boosting
        auto do_boost = [] (const TLorentzVector& bachelor, const TLorentzVector& etaprime) {
            TLorentzVector boosted(bachelor);
            boosted.Boost(-etaprime.BoostVector());
            return boosted;
        };

        t.Bachelor_E = do_boost(*g_EtaPrime_best, EtaPrime).E();

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

EtapOmegaG::Ref_t::Ref_t()
{

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
