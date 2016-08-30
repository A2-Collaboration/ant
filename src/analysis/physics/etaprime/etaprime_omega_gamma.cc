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
    params(utils::UncertaintyModels::Interpolated::makeAndLoad(
               // use FitterSergey as default
               make_shared<utils::UncertaintyModels::FitterSergey>(),
               utils::UncertaintyModels::Interpolated::Mode_t::Fit
               ),
           true, // flag to enable z vertex
           3.0 // Z_vertex_sigma, =0 means unmeasured
           ),
    kinfitter_sig("kinfitter_sig",4,
                  params.Fit_uncertainty_model, params.Fit_Z_vertex,
                  EtapOmegaG::MakeFitSettings(15)
                  ),
    kinfitter_ref("kinfitter_ref",2,
                  params.Fit_uncertainty_model, params.Fit_Z_vertex,
                  EtapOmegaG::MakeFitSettings(15)
                  ),
    mc_smear(opts->Get<bool>("MCFake", false) | opts->Get<bool>("MCSmear", true)
             ? // use | to force evaluation of both opts!
               std_ext::make_unique<utils::MCSmear>(
                   opts->Get<bool>("MCFake", false) ?
                       params.Fit_uncertainty_model // in Fake mode use same model as fitter
                     : utils::UncertaintyModels::Interpolated::makeAndLoad(
                           // use Adlarson as default (30% version of Oli is maybe better?)
                           make_shared<utils::UncertaintyModels::MCSmearingAdlarson>(),
                           utils::UncertaintyModels::Interpolated::Mode_t::MCSmear
                           )
                       )
             : nullptr // no MCSmear
               ),
    mc_fake(opts->Get<bool>("MCFake", false) ?
                std_ext::make_unique<utils::MCFakeReconstructed>()
              : nullptr),
    Sig(params)
{
    if(mc_smear)
        LOG(INFO) << "Additional MC Smearing enabled";

    const interval<double> prompt_range{-3.0,1.5};
    promptrandom.AddPromptRange(prompt_range); // slight offset due to CBAvgTime reference
    promptrandom.AddRandomRange({-35,-10});  // just ensure to be way off prompt peak
    promptrandom.AddRandomRange({ 10, 35});

    h_CommonCuts = HistFac.makeTH1D("Common Cuts", "", "#", BinSettings(15),"h_CommonCuts");
    h_CommonCuts_sig = HistFac.makeTH1D("Common Cuts Sig", "", "#", BinSettings(15),"h_CommonCuts_sig");
    h_CommonCuts_ref = HistFac.makeTH1D("Common Cuts Ref", "", "#", BinSettings(15),"h_CommonCuts_ref");
    h_MissedBkg = HistFac.makeTH1D("Missed Background", "", "#", BinSettings(25),"h_MissedBkg");

    h_LostPhotons_sig = HistFac.makeTH1D("LostPhotons Sig", "#theta", "#", BinSettings(200,0,180),"h_LostPhotons_sig");
    h_LostPhotons_ref = HistFac.makeTH1D("LostPhotons Ref", "#theta", "#", BinSettings(200,0,180),"h_LostPhotons_ref");

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

    const bool have_MCTrue = !event.MCTrue().ID.IsInvalid();

    const TEventData& data = mc_fake && have_MCTrue ? mc_fake->Get(event.MCTrue()) : event.Reconstructed();

    const bool is_MC = data.ID.isSet(TID::Flags_t::MC);

    h_CommonCuts->Fill("Seen",1.0);

    auto& particletree = event.MCTrue().ParticleTree;

    // Count EtaPrimes in MC sample
    h_CommonCuts->Fill("MCTrue #eta'", 0); // ensure the bin is there...
    if(particletree) {
        // note: this might also match to g p -> eta' eta' p,
        // but this is kinematically forbidden
        if(utils::ParticleTools::FindParticle(ParticleTypeDatabase::EtaPrime, particletree, 1)) {
            h_CommonCuts->Fill("MCTrue #eta'", 1);
        }
    }

    // do some MCTrue identification (if available)
    t.MCTrue = 0; // indicate data by default
    t.TrueZVertex = event.MCTrue().Target.Vertex.z; // NaN in case of data
    TParticleTree_t ptree_sig = nullptr; // used by Sig_t::Process to check matching
    if(particletree) {
        // 1=Signal, 2=Reference, 9=MissedBkg, >=10 found in ptreeBackgrounds
        if(particletree->IsEqual(ptreeSignal, utils::ParticleTools::MatchByParticleName)) {
            t.MCTrue = 1;
            ptree_sig = particletree;
        }
        else if(particletree->IsEqual(ptreeReference, utils::ParticleTools::MatchByParticleName)) {
            t.MCTrue = 2;
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
    else if(have_MCTrue) {
        // in rare cases, the particletree is not available, although we're running on MCTrue
        // mark this as other MC background
        t.MCTrue = 9;
    }

    // do some additional counting if true signal/ref event
    if(t.MCTrue == 1 || t.MCTrue == 2) {
        auto h_cut = t.MCTrue == 1 ? h_CommonCuts_sig : h_CommonCuts_ref;
        auto h_lost = t.MCTrue == 1 ? h_LostPhotons_sig : h_LostPhotons_ref;
        h_cut->Fill("MCTrue seen", 1.0);
        bool photons_accepted = true;
        for(const TParticlePtr& p : event.MCTrue().Particles.Get(ParticleTypeDatabase::Photon)) {
            if(geometry.DetectorFromAngles(*p) == Detector_t::Any_t::None) {
                h_lost->Fill(std_ext::radian_to_degree(p->Theta()));
                photons_accepted = false;
            }
        }
        if(photons_accepted) {
            h_cut->Fill("MCTrue Photon ok", 1.0);
        }
        auto proton = event.MCTrue().Particles.Get(ParticleTypeDatabase::Photon).front();
        if(geometry.DetectorFromAngles(*proton) != Detector_t::Any_t::None)
            h_cut->Fill("MCTrue Proton ok", 1.0);
    }

    // start now with some cuts

    // very simple trigger simulation for MC
    /// \todo Investigate trigger behaviour with pi0pi0 sample?
    if(is_MC) {
        if(data.Trigger.CBEnergySum<=550)
            return;
        h_CommonCuts->Fill("MC CBEnergySum>550",1.0);
    }
    t.CBSumE = data.Trigger.CBEnergySum;

    t.CBAvgTime = data.Trigger.CBTiming;
    if(!isfinite(t.CBAvgTime))
        return;
    h_CommonCuts->Fill("CBAvgTime ok",1.0);

    if(data.Candidates.size()<3)
        return;
    h_CommonCuts->Fill("nCands>=3", 1.0);

    // gather candidates sorted by energy
    TCandidatePtrList candidates;
    TCandidatePtrList candidates_taps;
    for(const auto& cand : data.Candidates.get_iter()) {
        if(cand->Detector & Detector_t::Type_t::TAPS) {
            candidates_taps.emplace_back(cand);
        }
        candidates.emplace_back(cand);
    }
    if(candidates_taps.empty())
        return;
    h_CommonCuts->Fill("1 in TAPS",1.0);

    std::sort(candidates.begin(), candidates.end(),
              [] (const TCandidatePtr& a, const TCandidatePtr& b) {
        return a->CaloEnergy > b->CaloEnergy;
    });

    // sum up the PID energy
    // (might be different to matched CB/PID Veto energy)
    t.PIDSumE = 0;
    for(const TCluster& cl : data.Clusters) {
        if(cl.DetectorType == Detector_t::Type_t::PID) {
            t.PIDSumE += cl.Energy;
        }
    }

    t.TaggNPrompt = 0;
    t.TaggNRandom = 0;

    auto true_proton = utils::ParticleTools::FindParticle(ParticleTypeDatabase::Proton, particletree);

    for(const TTaggerHit& taggerhit : data.TaggerHits) {
        promptrandom.SetTaggerHit(taggerhit.Time - t.CBAvgTime);
        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        if(promptrandom.State() == PromptRandom::Case::Prompt)
            t.TaggNPrompt()++;
        if(promptrandom.State() == PromptRandom::Case::Random)
            t.TaggNRandom()++;

        t.TaggW =  promptrandom.FillWeight();
        t.TaggE =  taggerhit.PhotonEnergy;
        t.TaggT =  taggerhit.Time;
        t.TaggCh = taggerhit.Channel;


        Sig.t.KinFitProb = std_ext::NaN;
        Ref.t.KinFitProb = std_ext::NaN;

        Particles_t best_sig_particles;
        Particles_t best_ref_particles;

        for(const auto& cand_proton :  candidates_taps) {

            Particles_t sig_particles;
            Particles_t ref_particles;

            auto proton = make_shared<TParticle>(ParticleTypeDatabase::Proton, cand_proton);
            sig_particles.Proton = proton;
            ref_particles.Proton = proton;

            for(const auto& cand : candidates) {
                if(cand == cand_proton)
                    continue;
                auto photon = make_shared<TParticle>(ParticleTypeDatabase::Photon, cand);

                if(sig_particles.Photons.size()<4) {
                    sig_particles.PhotonSum += *photon;
                    sig_particles.Photons.emplace_back(photon);
                }
                else
                    sig_particles.DiscardedEk += cand->CaloEnergy;

                if(ref_particles.Photons.size()<2) {
                    ref_particles.PhotonSum += *photon;
                    ref_particles.Photons.emplace_back(photon);
                }
                else
                    ref_particles.DiscardedEk += cand->CaloEnergy;
            }

            bool haveSig = sig_particles.Photons.size() >= 4;
            bool haveRef = ref_particles.Photons.size() >= 2;

            // additionally smear the particles in MC
            if(mc_smear && is_MC) {
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

            if(haveSig && doKinfit(taggerhit, true_proton, kinfitter_sig, sig_particles, Sig.t, h_CommonCuts_sig)) {
                best_sig_particles = sig_particles;
            }

            if(haveRef && doKinfit(taggerhit, true_proton, kinfitter_ref, ref_particles, Ref.t, h_CommonCuts_ref)) {
                best_ref_particles = ref_particles;
            }
        }

        Sig.ResetBranches();
        Ref.ResetBranches();

        bool dofill = false;

        if(!best_sig_particles.Photons.empty()) {
            Sig.Process(best_sig_particles, ptree_sig);
            dofill = true;
        }

        if(!best_ref_particles.Photons.empty()) {
            Ref.Process(best_ref_particles);
            dofill = true;
        }

        if(dofill) {
            t.Tree->Fill();
            Sig.Fill();
            Ref.t.Tree->Fill();
        }
    }

}

bool EtapOmegaG::doKinfit(const TTaggerHit& taggerhit,
                          TParticlePtr true_proton,
                          utils::KinFitter& kinfitter,
                          EtapOmegaG::Particles_t& particles,
                          EtapOmegaG::SharedTree_t& t,
                          TH1D* h_CommonCuts)
{
    h_CommonCuts->Fill("Seen KinFit", 1.0);

    const LorentzVec& photon_sum = particles.PhotonSum;

    // missing mass
    const LorentzVec beam_target = taggerhit.GetPhotonBeam() + LorentzVec({0, 0, 0}, ParticleTypeDatabase::Proton.Mass());
    const auto& missing_mass = (beam_target - photon_sum).M();

    const auto& missingmass_cut = ParticleTypeDatabase::Proton.GetWindow(400);
    if(!missingmass_cut.Contains(missing_mass))
        return false;
    h_CommonCuts->Fill("MM ok", 1.0);

    if(photon_sum.M()<600)
        return false;
    h_CommonCuts->Fill("IM ok", 1.0);

    if(particles.DiscardedEk>70)
        return false;
    h_CommonCuts->Fill("DiscEk ok", 1.0);

    particles.PhotonEnergy = taggerhit.PhotonEnergy;
    kinfitter.SetEgammaBeam(particles.PhotonEnergy);
    kinfitter.SetProton(particles.Proton);
    kinfitter.SetPhotons(particles.Photons);

    auto result = kinfitter.DoFit();

    if(result.Status != APLCON::Result_Status_t::Success)
        return false;

    if(!std_ext::copy_if_greater(t.KinFitProb, result.Probability))
        return false;

    if(t.KinFitProb<0.01)
        return false;

    t.DiscardedEk = particles.DiscardedEk;

    TParticlePtr& proton = particles.Proton;

    t.ProtonTime = proton->Candidate->Time;
    t.ProtonE = proton->Ek();
    t.ProtonTheta = std_ext::radian_to_degree(proton->Theta());
    t.ProtonVetoE = proton->Candidate->VetoEnergy;
    t.ProtonShortE = proton->Candidate->FindCaloCluster()->ShortEnergy;

    t.PhotonsEk = 0;
    t.nPhotonsCB = 0;
    t.nPhotonsTAPS = 0;
    t.CBSumVetoE = 0;
    t.PhotonThetas().clear();
    for(const auto& photon : particles.Photons) {
        const auto& cand = photon->Candidate;
        t.PhotonsEk += cand->CaloEnergy;
        if(cand->Detector & Detector_t::Type_t::CB) {
            t.nPhotonsCB++;
            t.CBSumVetoE += cand->VetoEnergy;
        }
        if(cand->Detector & Detector_t::Type_t::TAPS)
            t.nPhotonsTAPS++;
        t.PhotonThetas().emplace_back(std_ext::radian_to_degree(cand->Theta));
    }
    t.PhotonSum = photon_sum.M();

    t.ProtonCopl = std_ext::radian_to_degree(vec2::Phi_mpi_pi(proton->Phi() - photon_sum.Phi() - M_PI ));

    if(true_proton)
        t.ProtonTrueAngle = std_ext::radian_to_degree(proton->Angle(*true_proton));

    t.MissingMass = missing_mass;

    t.KinFitProb = result.Probability;
    t.KinFitIterations = result.NIterations;
    t.KinFitZVertex = kinfitter.GetFittedZVertex();

    t.KinFitBeamEPull = kinfitter.GetBeamEPull();
    t.KinFitProtonPulls = kinfitter.GetProtonPulls();
    t.KinFitPhotonsPulls = kinfitter.GetPhotonsPulls();

    const auto& fitted_proton = kinfitter.GetFittedProton();
    t.FittedProtonE = fitted_proton->Ek();

    particles.FittedPhotons = kinfitter.GetFittedPhotons();
    particles.FittedPhotonSum = {{0,0,0},0};

    for(const auto& photon : particles.FittedPhotons)
        particles.FittedPhotonSum += *photon;

    h_CommonCuts->Fill("KinFit ok", 1.0);

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
    AntiPi0ProtonPulls().resize(0);
    AntiPi0PhotonsPulls().resize(0);


    AntiEtaFitProb = std_ext::NaN;
    AntiEtaFitIterations = 0;
    AntiEtaFitZVertex = std_ext::NaN;

    AntiEtaBeamEPull = std_ext::NaN;
    AntiEtaProtonPulls().resize(0);
    AntiEtaPhotonsPulls().resize(0);
}

void EtapOmegaG::Sig_t::Process(const Particles_t& particles, const TParticleTree_t& ptree_sig)
{
    DoPhotonCombinatorics(particles.FittedPhotons);
    DoAntiPi0Eta(particles);
    if(t.AntiEtaFitProb>0.1 || t.AntiPi0FitProb>0.1)
        return;
    OmegaPi0.Process(particles, ptree_sig);
    Pi0.Process(particles, ptree_sig);
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
        t.AntiPi0ProtonPulls = fitter.GetProtonPulls();
        t.AntiPi0PhotonsPulls = fitter.GetPhotonsPulls();
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
        t.AntiEtaProtonPulls = fitter.GetProtonPulls();
        t.AntiEtaPhotonsPulls = fitter.GetPhotonsPulls();
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
    TreeFitProtonPulls().resize(0);
    TreeFitPhotonsPulls().resize(0);

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
                                       const TParticleTree_t& ptree_sig)
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
        t.TreeFitProtonPulls = treefitter.GetProtonPulls();
        t.TreeFitPhotonsPulls = treefitter.GetPhotonsPulls();

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
        if(ptree_sig) {
            auto true_photons = utils::ParticleTools::FindParticles(ParticleTypeDatabase::Photon, ptree_sig);
            assert(true_photons.size() == 4);
            auto match_bycandidate = [] (const TParticlePtr& mctrue, const TParticlePtr& recon) {
                return mctrue->Angle(*recon->Candidate); // TCandidate converts into vec3
            };
            auto matched = utils::match1to1(true_photons, particles.FittedPhotons,
                                            match_bycandidate,IntervalD(0.0, std_ext::degree_to_radian(15.0)));
            if(matched.size() == 4) {
                // find the two photons of the pi0
                TParticleList pi0_photons;
                ptree_sig->Map_nodes([&pi0_photons] (const TParticleTree_t& t) {
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
                                            const TParticleTree_t& ptree_sig)
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
        t.TreeFitProtonPulls = treefitter.GetProtonPulls();
        t.TreeFitPhotonsPulls = treefitter.GetPhotonsPulls();

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
        if(ptree_sig) {

            auto true_photons = utils::ParticleTools::FindParticles(ParticleTypeDatabase::Photon, ptree_sig);
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

                auto etap = select_daughter(ptree_sig, ParticleTypeDatabase::EtaPrime);
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
    canvas("Overview") << h_CommonCuts << h_MissedBkg
                       << h_CommonCuts_sig << h_CommonCuts_ref
                       << h_LostPhotons_sig << h_LostPhotons_ref
                       << endc;
    Ref.t.Tree->AddFriend(t.Tree);
    Sig.Pi0.t.Tree->AddFriend(t.Tree);
    Sig.Pi0.t.Tree->AddFriend(Sig.t.Tree);
    Sig.OmegaPi0.t.Tree->AddFriend(t.Tree);
    Sig.OmegaPi0.t.Tree->AddFriend(Sig.t.Tree);
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
