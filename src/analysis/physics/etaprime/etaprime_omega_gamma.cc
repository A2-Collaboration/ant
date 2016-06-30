#include "etaprime_omega_gamma.h"

#include "plot/root_draw.h"

#include "utils/particle_tools.h"
#include "utils/matcher.h"
#include "utils/combinatorics.h"
#include "analysis/utils/MCFakeReconstructed.h"

#include "expconfig/ExpConfig.h"

#include "base/Logger.h"
#include "base/std_ext/math.h"
#include "base/std_ext/vector.h"
#include "base/std_ext/misc.h"

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
//    settings.ConstraintAccuracy = 1.0e-3;
//    settings.Chi2Accuracy = 1.0e-2;
//    settings.DebugLevel = 5;
    return settings;
}

EtapOmegaG::EtapOmegaG(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    fit_Z_vertex(opts->Get<bool>("FitZVertex", true)),
    params(fit_Z_vertex ? make_shared<utils::UncertaintyModels::Optimized_Andi1>() :
                          make_shared<utils::UncertaintyModels::Optimized_Oli1>(),
           fit_Z_vertex, // flag to enable z vertex
           0.0 // Z_vertex_sigma, =0 means unmeasured
           ),
    kinfitter_sig("kinfitter_sig",4,
                  params.Fit_uncertainty_model, params.Fit_Z_vertex,
                  EtapOmegaG::MakeFitSettings(25)
                  ),
    kinfitter_ref("kinfitter_ref",2,
                  params.Fit_uncertainty_model, params.Fit_Z_vertex,
                  EtapOmegaG::MakeFitSettings(25)
                  ),
    mc_smear(opts->Get<bool>("MCFake", false) | opts->Get<bool>("MCSmear", true) ? // use | to force evaluation of both opts!
                 std_ext::make_unique<utils::MCSmear>(
                     opts->Get<bool>("MCFake", false) ?
                     params.Fit_uncertainty_model
                     : make_shared<utils::UncertaintyModels::MCSmearingAdlarson>()
                         )
               : nullptr
                 ),
    mc_fake(opts->Get<bool>("MCFake", false) ?
                std_ext::make_unique<utils::MCFakeReconstructed>()
              : nullptr),
    Sig(params)
{
    if(mc_smear)
        LOG(INFO) << "Additional MC Smearing enabled";

    const interval<double> prompt_range{-2.5,1.5};
    promptrandom.AddPromptRange(prompt_range); // slight offset due to CBAvgTime reference
    promptrandom.AddRandomRange({-30,-10});  // just ensure to be way off prompt peak
    promptrandom.AddRandomRange({ 10, 30});

    promptrandom_tight.AddPromptRange(prompt_range); // slight offset due to CBAvgTime reference
    promptrandom_tight.AddRandomRange({-20,-10});  // just ensure to be way off prompt peak
    promptrandom_tight.AddRandomRange({ 10, 20});


    h_CommonCuts = HistFac.makeTH1D("Common Cuts", "", "#", BinSettings(15),"h_CommonCuts");
    h_CommonCuts_sig = HistFac.makeTH1D("Common Cuts Sig", "", "#", BinSettings(15),"h_CommonCuts_sig");
    h_CommonCuts_ref = HistFac.makeTH1D("Common Cuts Ref", "", "#", BinSettings(15),"h_CommonCuts_ref");
    h_MissedBkg = HistFac.makeTH1D("Missed Background", "", "#", BinSettings(25),"h_MissedBkg");

    t.CreateBranches(HistFac.makeTTree("treeCommon"));

    Sig.SetupTrees(HistFac);
    Ref.t.CreateBranches(HistFac.makeTTree("Ref"));

    if(params.Fit_Z_vertex) {
        LOG(INFO) << "Fit Z vertex enabled with sigma=" << params.Z_vertex_sigma;
        kinfitter_sig.SetZVertexSigma(params.Z_vertex_sigma);
        kinfitter_ref.SetZVertexSigma(params.Z_vertex_sigma);
    }
}

void EtapOmegaG::ProcessEvent(const TEvent& event, manager_t&)
{
    // we start with some general candidate handling,
    // later we split into ref/sig analysis according to
    // number of photons

    const TEventData& data = mc_fake && !event.MCTrue().ID.IsInvalid() ? mc_fake->Get(event.MCTrue()) : event.Reconstructed();

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

    t.CBAvgTime = data.Trigger.CBTiming;
    if(!isfinite(t.CBAvgTime))
        return;
    h_CommonCuts->Fill("CBAvgTime ok",1.0);

    // gather candidates sorted by energy
    TCandidatePtrList candidates;
    bool haveTAPS = false;
    for(const auto& cand : data.Candidates.get_iter()) {
        if(cand->Detector & Detector_t::Type_t::TAPS) {
            haveTAPS = true;
        }
        candidates.emplace_back(cand);
    }
    if(!haveTAPS)
        return;
    h_CommonCuts->Fill("1 in TAPS",1.0);

    std::sort(candidates.begin(), candidates.end(),
              [] (const TCandidatePtr& a, const TCandidatePtr& b) {
        return a->CaloEnergy > b->CaloEnergy;
    });

    // fill sig/ref particles
    Particles_t sig_particles;
    Particles_t ref_particles;
    auto true_proton = utils::ParticleTools::FindParticle(ParticleTypeDatabase::Proton, particletree);
    bool haveSig = findParticles(candidates, 4, true_proton, sig_particles, Sig.t, h_CommonCuts_sig);
    bool haveRef = findParticles(candidates, 2, true_proton, ref_particles, Ref.t, h_CommonCuts_ref);

    if(!haveSig && !haveRef)
        return;
    h_CommonCuts->Fill("Sig||Ref",1.0);

    // sum up the PID energy
    // (might be different to matched CB/PID Veto energy)
    t.PIDSumE = 0;
    for(const TCluster& cl : data.Clusters) {
        if(cl.DetectorType == Detector_t::Type_t::PID) {
            t.PIDSumE += cl.Energy;
        }
    }

    // do some MCTrue identification (if available)
    t.MCTrue = 0; // indicate data by default
    t.TrueZVertex = event.MCTrue().Target.Vertex.z; // NaN in case of data
    TParticleTree_t ptree_sigref = nullptr; // used by Sig_t::Process to check matching
    if(particletree) {
        // 1=Signal, 2=Reference, 9=MissedBkg, >=10 found in ptreeBackgrounds
        if(particletree->IsEqual(ptreeSignal, utils::ParticleTools::MatchByParticleName)) {
            t.MCTrue = 1;
            ptree_sigref = particletree;
        }
        else if(particletree->IsEqual(ptreeReference, utils::ParticleTools::MatchByParticleName)) {
            t.MCTrue = 2;
            ptree_sigref = particletree;
        }
        else {
            t.MCTrue = 10;
            bool found = false;
            for(const auto& ptreeBkg : ptreeBackgrounds) {
                if(particletree->IsEqual(ptreeBkg.Tree, utils::ParticleTools::MatchByParticleName)) {
                    found = true;
                    break;
                }
                t.MCTrue++;
            }
            if(!found) {
                t.MCTrue = 9;
                const auto& decaystr = utils::ParticleTools::GetDecayString(particletree);
                h_MissedBkg->Fill(decaystr.c_str(), 1.0);
            }
        }
    }
    else if(!event.MCTrue().ID.IsInvalid()) {
        // in rare cases, the particletree is not available, although we're running on MCTrue
        // mark this as other MC background
        t.MCTrue = 9;
    }

    // additionally smear the particles in MC
    if(mc_smear && data.ID.isSet(TID::Flags_t::MC)) {
        auto smear_particles = [this] (Particles_t& particles) {
            particles.Proton = mc_smear->Smear(particles.Proton);
            for(auto& p : particles.Photons)
                p = mc_smear->Smear(p);
        };
        if(haveSig)
            smear_particles(sig_particles);
        if(haveRef)
            smear_particles(ref_particles);
    }

    for(const TTaggerHit& taggerhit : data.TaggerHits) {
        promptrandom.SetTaggerHit(taggerhit.Time - t.CBAvgTime);
        promptrandom_tight.SetTaggerHit(taggerhit.Time - t.CBAvgTime);
        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        t.TaggW = promptrandom.FillWeight();
        t.TaggW_tight = promptrandom_tight.FillWeight();
        t.TaggE = taggerhit.PhotonEnergy;
        t.TaggT = taggerhit.Time;
        t.TaggCh = taggerhit.Channel;

        Sig.ResetBranches();
        Ref.ResetBranches();

        bool dofill = false;

        if(haveSig && doKinfit(taggerhit, kinfitter_sig, sig_particles, Sig.t)) {
            Sig.Process(sig_particles, ptree_sigref);
            dofill = true;
        }

        if(haveRef && doKinfit(taggerhit, kinfitter_ref, ref_particles, Ref.t)) {
            Ref.Process(ref_particles);
            dofill = true;
        }

        if(dofill) {
            t.Tree->Fill();
            Sig.Fill();
            Ref.t.Tree->Fill();
        }
    }

}

bool EtapOmegaG::findParticles(const TCandidatePtrList& candidates,
                               unsigned nPhotons,
                               TParticlePtr true_proton,
                               EtapOmegaG::Particles_t& particles,
                               EtapOmegaG::SharedTree_t& t,
                               TH1D* h_CommonCuts)
{
    t.ProtonE = std_ext::NaN;
    t.ProtonTheta = std_ext::NaN;
    t.ProtonVetoE = std_ext::NaN;
    t.ProtonShortE = std_ext::NaN;
    t.ProtonTrueAngle = std_ext::NaN;
    t.DiscardedEk = std_ext::NaN;
    t.PhotonsEk = std_ext::NaN;
    t.nPhotonsCB   = 0;
    t.nPhotonsTAPS = 0;
    t.CBSumVetoE = std_ext::NaN;
    t.PhotonSum  = std_ext::NaN;
    t.PhotonThetas().resize(0);
    t.ProtonCopl = std_ext::NaN;

    t.nCandidates = candidates.size();
    if(t.nCandidates<nPhotons+1)
        return false;

    // identify the proton here as slowest cluster in TAPS
    // using CBAvgTime does not help here, since it's constant
    // over each TAPS clusters
    /// \todo think about using beta here as in EtapProton?
    t.ProtonTime = std_ext::NaN;
    TParticlePtr& proton = particles.Proton;
    for(unsigned i=0;i<nPhotons+1;i++) {
        const auto& cand = candidates[i];
        if(cand->Detector & Detector_t::Type_t::TAPS) {
            if(!isfinite(t.ProtonTime) || t.ProtonTime < cand->Time) {
                t.ProtonTime = cand->Time;
                proton = make_shared<TParticle>(ParticleTypeDatabase::Proton, cand);
            }
        }
    }

    if(!proton)
        return false;
    h_CommonCuts->Fill("p in TAPS", 1.0);

    t.ProtonE = proton->Ek();
    t.ProtonTheta = std_ext::radian_to_degree(proton->Theta());
    t.ProtonVetoE = proton->Candidate->VetoEnergy;
    t.ProtonShortE = proton->Candidate->FindCaloCluster()->ShortEnergy;

    // sum up discarded candidates
    t.DiscardedEk = 0;
    for(unsigned i=nPhotons+1;i<t.nCandidates;i++) {
        const auto& cand = candidates[i];
        t.DiscardedEk += cand->CaloEnergy;
    }

    // remaining candidates are photons
    LorentzVec& photon_sum = particles.PhotonSum;
    photon_sum = {0,0,0,0};
    t.PhotonsEk = 0;
    t.nPhotonsCB = 0;
    t.nPhotonsTAPS = 0;
    t.CBSumVetoE = 0;
    for(unsigned i=0;i<nPhotons+1;i++) {
        const auto& cand = candidates[i];
        if(cand == proton->Candidate)
            continue;
        t.PhotonsEk += cand->CaloEnergy;
        if(cand->Detector & Detector_t::Type_t::CB) {
            t.nPhotonsCB++;
            t.CBSumVetoE += cand->VetoEnergy;
        }
        if(cand->Detector & Detector_t::Type_t::TAPS)
            t.nPhotonsTAPS++;
        auto photon = make_shared<TParticle>(ParticleTypeDatabase::Photon, cand);
        photon_sum += *photon;
        t.PhotonThetas().emplace_back(std_ext::radian_to_degree(cand->Theta));
        particles.Photons.emplace_back(move(photon));
    }
    assert(particles.Photons.size() == nPhotons);
    assert(nPhotons == t.nPhotonsCB + t.nPhotonsTAPS);

    t.PhotonSum = photon_sum.M();

    // don't bother with events where proton coplanarity is not ok
    // we use some rather large window here...
    t.ProtonCopl = std_ext::radian_to_degree(vec2::Phi_mpi_pi(proton->Phi() - photon_sum.Phi() - M_PI ));
    const interval<double> ProtonCopl_cut(-35, 35);
    if(!ProtonCopl_cut.Contains(t.ProtonCopl))
        return false;
    h_CommonCuts->Fill("ProtonCopl ok", 1.0);

    if(true_proton)
        t.ProtonTrueAngle = std_ext::radian_to_degree(proton->Angle(*true_proton));

    return true;
}

bool EtapOmegaG::doKinfit(const TTaggerHit& taggerhit,
                          utils::KinFitter& kinfitter,
                          EtapOmegaG::Particles_t& particles,
                          EtapOmegaG::SharedTree_t& t)
{
    // missing mass
    const LorentzVec beam_target = taggerhit.GetPhotonBeam() + LorentzVec(0, 0, 0, ParticleTypeDatabase::Proton.Mass());
    t.MissingMass = (beam_target - particles.PhotonSum).M();

    t.KinFitProb = std_ext::NaN;
    t.KinFitIterations = 0;
    t.KinFitZVertex = std_ext::NaN;

    t.KinFitBeamEPull = std_ext::NaN;
    t.KinFitProtonEPull = std_ext::NaN;
    t.KinFitProtonThetaPull = std_ext::NaN;
    t.KinFitProtonPhiPull = std_ext::NaN;
    t.KinFitPhotonEPulls().resize(0);
    t.KinFitPhotonThetaPulls().resize(0);
    t.KinFitPhotonPhiPulls().resize(0);

    t.FittedProtonE = std_ext::NaN;

    const auto& missingmass_cut = ParticleTypeDatabase::Proton.GetWindow(300);
    if(!missingmass_cut.Contains(t.MissingMass))
        return false;

    particles.PhotonEnergy = taggerhit.PhotonEnergy;
    kinfitter.SetEgammaBeam(particles.PhotonEnergy);
    kinfitter.SetProton(particles.Proton);
    kinfitter.SetPhotons(particles.Photons);

    auto result = kinfitter.DoFit();

    if(result.Status != APLCON::Result_Status_t::Success)
        return false;

    t.KinFitProb = result.Probability;
    t.KinFitIterations = result.NIterations;
    t.KinFitZVertex = kinfitter.GetFittedZVertex();

    t.KinFitBeamEPull = kinfitter.GetBeamEPull();
    t.KinFitProtonEPull = kinfitter.GetProtonEPull();
    t.KinFitProtonThetaPull = kinfitter.GetProtonThetaPull();
    t.KinFitProtonPhiPull = kinfitter.GetProtonPhiPull();
    t.KinFitPhotonEPulls = kinfitter.GetPhotonEPulls();
    t.KinFitPhotonThetaPulls = kinfitter.GetPhotonThetaPulls();
    t.KinFitPhotonPhiPulls = kinfitter.GetPhotonPhiPulls();

    const auto& fitted_proton = kinfitter.GetFittedProton();
    t.FittedProtonE = fitted_proton->Ek();

    particles.FittedPhotons = kinfitter.GetFittedPhotons();
    particles.FittedPhotonSum = {0,0,0,0};

    for(const auto& photon : particles.FittedPhotons)
        particles.FittedPhotonSum += *photon;

    return true;
}

EtapOmegaG::Sig_t::Sig_t(params_t params) :
    Pi0(params),
    OmegaPi0(params),
    treefitter_Pi0Pi0("treefit_Pi0Pi0",
                      ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::TwoPi0_4g),
                      params.Fit_uncertainty_model, params.Fit_Z_vertex, {},
                      MakeFitSettings(20)
                      ),
    treefitter_Pi0Eta("treefit_Pi0Eta",
                      ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0Eta_4g),
                      params.Fit_uncertainty_model, params.Fit_Z_vertex, {},
                      MakeFitSettings(20)
                      )
{
    if(params.Fit_Z_vertex) {
        treefitter_Pi0Pi0.SetZVertexSigma(params.Z_vertex_sigma);
        treefitter_Pi0Eta.SetZVertexSigma(params.Z_vertex_sigma);
    }
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


template<typename T>
void fill_NaN(T& v) {
    std::fill(v.begin(), v.end(), std_ext::NaN);
}

void EtapOmegaG::Sig_t::SharedTree_t::Reset()
{
    fill_NaN(ggg());
    fill_NaN(gg_gg1());
    fill_NaN(gg_gg2());

    AntiPi0FitProb = std_ext::NaN;
    AntiPi0FitIterations = 0;
    AntiPi0FitZVertex = std_ext::NaN;

    AntiPi0BeamEPull = std_ext::NaN;
    AntiPi0ProtonEPull = std_ext::NaN;
    AntiPi0ProtonThetaPull = std_ext::NaN;
    AntiPi0ProtonPhiPull = std_ext::NaN;
    AntiPi0PhotonEPulls().resize(0);
    AntiPi0PhotonThetaPulls().resize(0);
    AntiPi0PhotonPhiPulls().resize(0);

    AntiEtaFitProb = std_ext::NaN;
    AntiEtaFitIterations = 0;
    AntiEtaFitZVertex = std_ext::NaN;

    AntiEtaBeamEPull = std_ext::NaN;
    AntiEtaProtonEPull = std_ext::NaN;
    AntiEtaProtonThetaPull = std_ext::NaN;
    AntiEtaProtonPhiPull = std_ext::NaN;
    AntiEtaPhotonEPulls().resize(0);
    AntiEtaPhotonThetaPulls().resize(0);
    AntiEtaPhotonPhiPulls().resize(0);
}

void EtapOmegaG::Sig_t::Process(const Particles_t& particles, const TParticleTree_t& ptree_sigref)
{
    DoPhotonCombinatorics(particles.FittedPhotons);
    DoAntiPi0Eta(particles);
    OmegaPi0.Process(particles, ptree_sigref);
    Pi0.Process(particles, ptree_sigref);
}

void EtapOmegaG::Sig_t::DoPhotonCombinatorics(const TParticleList& photons)
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

void EtapOmegaG::Sig_t::DoAntiPi0Eta(const Particles_t& particles)
{
    APLCON::Result_t r;

    treefitter_Pi0Pi0.SetEgammaBeam(particles.PhotonEnergy);
    treefitter_Pi0Pi0.SetProton(particles.Proton);
    treefitter_Pi0Pi0.SetPhotons(particles.Photons);
    while(treefitter_Pi0Pi0.NextFit(r)) {
        if(r.Status != APLCON::Result_Status_t::Success)
            continue;
        if(!std_ext::copy_if_greater(t.AntiPi0FitProb, r.Probability))
            continue;
        // found fit with better prob
        t.AntiPi0FitIterations = r.NIterations;

        const auto& fitter = treefitter_Pi0Pi0;
        t.AntiPi0FitZVertex = fitter.GetFittedZVertex();
        t.AntiPi0BeamEPull = fitter.GetBeamEPull();
        t.AntiPi0ProtonEPull = fitter.GetProtonEPull();
        t.AntiPi0ProtonThetaPull = fitter.GetProtonThetaPull();
        t.AntiPi0ProtonPhiPull = fitter.GetProtonPhiPull();
        t.AntiPi0PhotonEPulls = fitter.GetPhotonEPulls();
        t.AntiPi0PhotonThetaPulls = fitter.GetPhotonThetaPulls();
        t.AntiPi0PhotonPhiPulls = fitter.GetPhotonPhiPulls();
    }

    treefitter_Pi0Eta.SetEgammaBeam(particles.PhotonEnergy);
    treefitter_Pi0Eta.SetProton(particles.Proton);
    treefitter_Pi0Eta.SetPhotons(particles.Photons);
    while(treefitter_Pi0Eta.NextFit(r)) {
        if(r.Status != APLCON::Result_Status_t::Success)
            continue;
        if(!std_ext::copy_if_greater(t.AntiEtaFitProb, r.Probability))
            continue;
        // found fit with better probability
        t.AntiEtaFitIterations = r.NIterations;

        const auto& fitter = treefitter_Pi0Eta;
        t.AntiEtaFitZVertex = fitter.GetFittedZVertex();
        t.AntiEtaBeamEPull = fitter.GetBeamEPull();
        t.AntiEtaProtonEPull = fitter.GetProtonEPull();
        t.AntiEtaProtonThetaPull = fitter.GetProtonThetaPull();
        t.AntiEtaProtonPhiPull = fitter.GetProtonPhiPull();
        t.AntiEtaPhotonEPulls = fitter.GetPhotonEPulls();
        t.AntiEtaPhotonThetaPulls = fitter.GetPhotonThetaPulls();
        t.AntiEtaPhotonPhiPulls = fitter.GetPhotonPhiPulls();
    }
}

EtapOmegaG::Sig_t::Fit_t::Fit_t(utils::TreeFitter fitter) :
    treefitter(move(fitter)),
    fitted_Pi0(treefitter.GetTreeNode(ParticleTypeDatabase::Pi0)),
    fitted_Omega(treefitter.GetTreeNode(ParticleTypeDatabase::Omega)),
    fitted_EtaPrime(treefitter.GetTreeNode(ParticleTypeDatabase::EtaPrime))
{

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

    fitted_g_Omega = find_photons(fitted_Omega).at(0);

    fitted_g_EtaPrime = find_photons(fitted_EtaPrime).at(0);
}

utils::TreeFitter EtapOmegaG::Sig_t::Fit_t::Make(const ParticleTypeDatabase::Type& subtree, params_t params)
{
    auto setupnodes = [&subtree] (const ParticleTypeTree& t) {
        utils::TreeFitter::nodesetup_t nodesetup;
        // always exlude the EtaPrime
        if(t->Get() == ParticleTypeDatabase::EtaPrime)
            nodesetup.Excluded = true;
        // subtree decides if the Omega is excluded or not
        if(subtree == ParticleTypeDatabase::Pi0 &&
           t->Get() == ParticleTypeDatabase::Omega)
            nodesetup.Excluded = true;
        return nodesetup;
    };

    utils::TreeFitter treefitter{
        "sig_treefitter_"+subtree.Name(),
        EtapOmegaG::ptreeSignal,
        params.Fit_uncertainty_model,
        params.Fit_Z_vertex,
        setupnodes,
        MakeFitSettings(15)
    };
    if(params.Fit_Z_vertex)
        treefitter.SetZVertexSigma(params.Z_vertex_sigma);
    return treefitter;
}

TParticlePtr EtapOmegaG::Sig_t::Fit_t::FindBest(const utils::TreeFitter::tree_t& fitted,
                                                const Particles_t& particles)
{
    // use at() to have boundary check
    return particles.FittedPhotons.at(fitted->Get().PhotonLeaveIndex);
}

void EtapOmegaG::Sig_t::Fit_t::BaseTree_t::Reset()
{
    TreeFitProb = std_ext::NaN;
    TreeFitIterations = 0;
    TreeFitZVertex = std_ext::NaN;

    TreeFitBeamEPull = std_ext::NaN;
    TreeFitProtonEPull = std_ext::NaN;
    TreeFitProtonThetaPull = std_ext::NaN;
    TreeFitProtonPhiPull = std_ext::NaN;
    TreeFitPhotonEPulls().resize(0);
    TreeFitPhotonThetaPulls().resize(0);
    TreeFitPhotonPhiPulls().resize(0);

    IM_Pi0_fitted = std_ext::NaN;
    IM_Pi0_best = std_ext::NaN;
    IM_Pi0gg_fitted = std_ext::NaN;
    IM_Pi0gg_best = std_ext::NaN;
    IM_gg_fitted = std_ext::NaN;
    IM_gg_best = std_ext::NaN;


    MCTrueMatch = 0;
}

EtapOmegaG::Sig_t::Pi0_t::Pi0_t(params_t params) :
    Fit_t(Fit_t::Make(ParticleTypeDatabase::Pi0, params))
{

}

void EtapOmegaG::Sig_t::Pi0_t::BaseTree_t::Reset()
{
    Fit_t::BaseTree_t::Reset();

    fill_NaN(IM_Pi0g_fitted());
    fill_NaN(IM_Pi0g_best());

    fill_NaN(Bachelor_E_fitted());
    fill_NaN(Bachelor_E_best());

}

void EtapOmegaG::Sig_t::Pi0_t::Process(const EtapOmegaG::Particles_t& particles,
                                       const TParticleTree_t& ptree_sigref)
{
    assert(particles.Photons.size() == 4);

    // the EtaPrime
    t.IM_Pi0gg_best = particles.FittedPhotonSum.M();

    // for MCtrue identification
    TParticlePtr g1_Pi0_best;
    TParticlePtr g2_Pi0_best;

    // do treefit
    treefitter.SetEgammaBeam(particles.PhotonEnergy);
    treefitter.SetProton(particles.Proton);
    treefitter.SetPhotons(particles.Photons);

    APLCON::Result_t r;

    while(treefitter.NextFit(r)) {
        if(r.Status != APLCON::Result_Status_t::Success)
            continue;
        if(!std_ext::copy_if_greater(t.TreeFitProb, r.Probability))
            continue;
        // found fit with better prob
        t.TreeFitIterations = r.NIterations;

        t.TreeFitZVertex = treefitter.GetFittedZVertex();
        t.TreeFitBeamEPull = treefitter.GetBeamEPull();
        t.TreeFitProtonEPull = treefitter.GetProtonEPull();
        t.TreeFitProtonThetaPull = treefitter.GetProtonThetaPull();
        t.TreeFitProtonPhiPull = treefitter.GetProtonPhiPull();
        t.TreeFitPhotonEPulls = treefitter.GetPhotonEPulls();
        t.TreeFitPhotonThetaPulls = treefitter.GetPhotonThetaPulls();
        t.TreeFitPhotonPhiPulls = treefitter.GetPhotonPhiPulls();

        // IM fitted expected to be delta peaks since they were fitted...
        const LorentzVec& Pi0_fitted = fitted_Pi0->Get().LVSum;
        t.IM_Pi0_fitted = Pi0_fitted.M();

        t.IM_Pi0gg_fitted = fitted_EtaPrime->Get().LVSum.M();

        // have a look at the assigned gammas to Pi0
        g1_Pi0_best = FindBest(fitted_g1_Pi0, particles);
        g2_Pi0_best = FindBest(fitted_g2_Pi0, particles);

        const LorentzVec& Pi0_best = *g1_Pi0_best + *g2_Pi0_best;
        t.IM_Pi0_best = Pi0_best.M();

        // there are two photon combinations possible
        // for the omega
        auto do_comb = [] (
                       const TParticlePtr& g1, const TParticlePtr& g2,
                       const LorentzVec Pi0,
                       double& IM_gg,
                       vector<double>& IM_Pi0g,
                       vector<double>& Bachelor_E
                       )
        {
            assert(IM_Pi0g.size() == 2);
            IM_gg = (*g1 + *g2).M();
            IM_Pi0g.front() = (Pi0 + *g1).M();
            IM_Pi0g.back()  = (Pi0 + *g2).M();
            const LorentzVec EtaPrime = *g1 + *g2 + Pi0;
            Bachelor_E.front() = Boost(*g1, -EtaPrime.BoostVector()).E;
            Bachelor_E.back() =  Boost(*g2, -EtaPrime.BoostVector()).E;

            if(IM_Pi0g.front() > IM_Pi0g.back()) {
                std::swap(IM_Pi0g.front(), IM_Pi0g.back());
                std::swap(Bachelor_E.front(), Bachelor_E.back());
            }
        };

        do_comb(FindBest(fitted_g_Omega, particles),
                FindBest(fitted_g_EtaPrime, particles),
                Pi0_best,
                t.IM_gg_best,
                t.IM_Pi0g_best,
                t.Bachelor_E_best);


        do_comb(fitted_g_Omega->Get().Leave->AsFitted(),
                fitted_g_EtaPrime->Get().Leave->AsFitted(),
                Pi0_fitted,
                t.IM_gg_fitted,
                t.IM_Pi0g_fitted,
                t.Bachelor_E_fitted);

    }


    // there was at least one successful fit
    if(isfinite(t.TreeFitProb)) {

        // check MC matching
        if(ptree_sigref) {
            auto true_photons = utils::ParticleTools::FindParticles(ParticleTypeDatabase::Photon, ptree_sigref);
            assert(true_photons.size() == 4);
            auto match_bycandidate = [] (const TParticlePtr& mctrue, const TParticlePtr& recon) {
                return mctrue->Angle(*recon->Candidate); // TCandidate converts into vec3
            };
            auto matched = utils::match1to1(true_photons, particles.FittedPhotons,
                                            match_bycandidate,IntervalD(0.0, std_ext::degree_to_radian(15.0)));
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


EtapOmegaG::Sig_t::OmegaPi0_t::OmegaPi0_t(params_t params) :
    Fit_t(Fit_t::Make(ParticleTypeDatabase::Omega, params))
{

}

void EtapOmegaG::Sig_t::OmegaPi0_t::BaseTree_t::Reset()
{
    Fit_t::BaseTree_t::Reset();
    IM_Pi0g_fitted = std_ext::NaN;
    IM_Pi0g_best = std_ext::NaN;
    Bachelor_E_fitted = std_ext::NaN;
    Bachelor_E_best = std_ext::NaN;

}


void EtapOmegaG::Sig_t::OmegaPi0_t::Process(const EtapOmegaG::Particles_t& particles,
                                            const TParticleTree_t& ptree_sigref)
{

    assert(particles.Photons.size() == 4);


    // g_Omega to check against MCTrue
    TParticlePtr g_Omega_best;
    // the EtaPrime bachelor photon is most important to us...
    TParticlePtr g_EtaPrime_best;
    TParticlePtr g_EtaPrime_fitted;
    const LorentzVec& EtaPrime_best = particles.FittedPhotonSum;
    LorentzVec EtaPrime_fitted;
    t.IM_Pi0gg_best = EtaPrime_best.M();


    // do treefit
    treefitter.SetEgammaBeam(particles.PhotonEnergy);
    treefitter.SetProton(particles.Proton);
    treefitter.SetPhotons(particles.Photons);

    APLCON::Result_t r;

    while(treefitter.NextFit(r)) {
        if(r.Status != APLCON::Result_Status_t::Success)
            continue;
        if(!std_ext::copy_if_greater(t.TreeFitProb, r.Probability))
            continue;
        // found fit with better prob
        t.TreeFitIterations = r.NIterations;
        t.TreeFitZVertex = treefitter.GetFittedZVertex();

        t.TreeFitBeamEPull = treefitter.GetBeamEPull();
        t.TreeFitProtonEPull = treefitter.GetProtonEPull();
        t.TreeFitProtonThetaPull = treefitter.GetProtonThetaPull();
        t.TreeFitProtonPhiPull = treefitter.GetProtonPhiPull();
        t.TreeFitPhotonEPulls = treefitter.GetPhotonEPulls();
        t.TreeFitPhotonThetaPulls = treefitter.GetPhotonThetaPulls();
        t.TreeFitPhotonPhiPulls = treefitter.GetPhotonPhiPulls();

        // IM fitted expected to be delta peaks since they were fitted...
        EtaPrime_fitted = fitted_EtaPrime->Get().LVSum;
        t.IM_Pi0gg_fitted = EtaPrime_fitted.M();
        t.IM_Pi0g_fitted = fitted_Omega->Get().LVSum.M();
        t.IM_Pi0_fitted = fitted_Pi0->Get().LVSum.M();


        // have a look at the assigned gammas to Pi0/Omega
        const auto& g1_Pi0_best = FindBest(fitted_g1_Pi0, particles);
        const auto& g2_Pi0_best = FindBest(fitted_g2_Pi0, particles);

        const LorentzVec& Pi0_best = *g1_Pi0_best + *g2_Pi0_best;
        t.IM_Pi0_best = Pi0_best.M();

        g_Omega_best = FindBest(fitted_g_Omega, particles);
        const LorentzVec& Omega_best = *g_Omega_best + Pi0_best;
        t.IM_Pi0g_best = Omega_best.M();

        // have a look at the EtaPrime bachelor photon
        // the element NOT in the combination is the Bachelor photon
        g_EtaPrime_best = FindBest(fitted_g_EtaPrime, particles);
        g_EtaPrime_fitted = fitted_g_EtaPrime->Get().Leave->AsFitted();


        t.IM_gg_best = (*g_EtaPrime_best + *g_Omega_best).M();
        t.IM_gg_fitted = (  *fitted_g_EtaPrime->Get().Leave->AsFitted()
                          + *fitted_g_Omega->Get().Leave->AsFitted()).M();

    }



    // there was at least one successful fit
    if(isfinite(t.TreeFitProb)) {

        t.Bachelor_E_best   = Boost(*g_EtaPrime_best,   -EtaPrime_best.BoostVector()).E;
        t.Bachelor_E_fitted = Boost(*g_EtaPrime_fitted, -EtaPrime_fitted.BoostVector()).E;

        // check MC matching
        if(ptree_sigref) {

            auto true_photons = utils::ParticleTools::FindParticles(ParticleTypeDatabase::Photon, ptree_sigref);
            assert(true_photons.size() == 4);
            auto match_bycandidate = [] (const TParticlePtr& mctrue, const TParticlePtr& recon) {
                return mctrue->Angle(*recon->Candidate); // TCandidate converts into vec3
            };
            auto matched = utils::match1to1(true_photons, particles.FittedPhotons,
                                            match_bycandidate,
                                            IntervalD(0.0, std_ext::degree_to_radian(15.0)));
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

void EtapOmegaG::Ref_t::Tree_t::Reset()
{
    IM_2g = std_ext::NaN;
}

void EtapOmegaG::Ref_t::ResetBranches()
{
    t.Reset();
}

void EtapOmegaG::Ref_t::Process(const EtapOmegaG::Particles_t& particles) {
    assert(particles.FittedPhotons.size() == 2);
    t.IM_2g = particles.FittedPhotonSum.M();
}

void EtapOmegaG::ShowResult()
{
    canvas("Overview") << h_CommonCuts  << h_CommonCuts_sig
                       << h_CommonCuts_ref << h_MissedBkg << endc;
//    if(fit_Z_vertex) {
        Ref.t.Tree->AddFriend(t.Tree);
        Sig.Pi0.t.Tree->AddFriend(t.Tree);
        Sig.Pi0.t.Tree->AddFriend(Sig.t.Tree);
        Sig.OmegaPi0.t.Tree->AddFriend(t.Tree);
        Sig.OmegaPi0.t.Tree->AddFriend(Sig.t.Tree);
//        canvas("Z Vertex Sig")
//                << TTree_drawable(Sig.Pi0.t.Tree, "KinFitZVertex:TrueZVertex >> h3(100,-5,5,100,-5,5)","KinFitProb>0.01")
//                << TTree_drawable(Sig.Pi0.t.Tree, "IM_Pi0gg_fitted:TrueZVertex >> h4(100,-5,5,200,900,1000)","KinFitProb>0.01")
//                << TTree_drawable(Sig.OmegaPi0.t.Tree, "KinFitZVertex:TrueZVertex >> h5(100,-5,5,100,-5,5)","KinFitProb>0.01")
//                << TTree_drawable(Sig.OmegaPi0.t.Tree, "IM_Pi0gg_fitted:TrueZVertex >> h6(100,-5,5,200,900,1000)","KinFitProb>0.01")
//                << endc;
        canvas("Z Vertex Ref")
                << drawoption("colz")
                << TTree_drawable(Ref.t.Tree, "KinFitZVertex:TrueZVertex >> h1(100,-5,5,100,-5,5)","KinFitProb>0.01")
                << TTree_drawable(Ref.t.Tree, "IM_2g:TrueZVertex >> h2(100,-5,5,200,900,1000)","KinFitProb>0.01")
                << TTree_drawable(Ref.t.Tree, "KinFitPhotonThetaPulls:PhotonThetas >> h3(100,5,175,50,-3,3)","KinFitProb>0.01")
                << endc;
//    }
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
