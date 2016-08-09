#include "etaprime_ept.h"

#include "plot/root_draw.h"

#include "utils/particle_tools.h"
#include "utils/matcher.h"
#include "utils/combinatorics.h"
#include "analysis/utils/MCFakeReconstructed.h"

#include "expconfig/ExpConfig.h"
#include "expconfig/detectors/EPT.h"

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

APLCON::Fit_Settings_t EtapEPT::MakeFitSettings(unsigned max_iterations)
{
    auto settings = APLCON::Fit_Settings_t::Default;
    settings.MaxIterations = max_iterations;
//    settings.ConstraintAccuracy = 1.0e-3;
//    settings.Chi2Accuracy = 1.0e-2;
//    settings.DebugLevel = 5;
    return settings;
}

EtapEPT::EtapEPT(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    EPT(ExpConfig::Setup::GetDetector<expconfig::detector::EPT>()),
    noTaggChPerm(opts->Get<bool>("NoTaggChPerm", false)),
    kinfitter("kinfitter",2,
              utils::UncertaintyModels::Interpolated::makeAndLoad(
                  // use OptimizedOli1 as default
                  make_shared<utils::UncertaintyModels::Optimized_Oli1>(),
                  utils::UncertaintyModels::Interpolated::Mode_t::Fit
                  ),
              true,
              EtapEPT::MakeFitSettings(25)
              ),
    mc_smear(opts->Get<bool>("MCSmear", true) ?
                 std_ext::make_unique<utils::MCSmear>(
                     utils::UncertaintyModels::Interpolated::makeAndLoad(
                         // use Adlarson as default (30% version of Oli is maybe better?)
                         make_shared<utils::UncertaintyModels::MCSmearingAdlarson>(),
                         utils::UncertaintyModels::Interpolated::Mode_t::MCSmear
                         )
                     )
               : nullptr // no MCSmear
             )
{
    if(mc_smear)
        LOG(INFO) << "Additional MC Smearing enabled";

    promptrandom.AddPromptRange({-2.5,1.5}); // slight offset due to CBAvgTime reference
    promptrandom.AddRandomRange({-30,-10});  // just ensure to be way off prompt peak
    promptrandom.AddRandomRange({ 10, 30});

    promptrandom_wide.AddPromptRange({ -7,  7});
    promptrandom_wide.AddRandomRange({-65,-15});  // just ensure to be way off prompt peak
    promptrandom_wide.AddRandomRange({ 15, 65});

    h_Cuts = HistFac.makeTH1D("Cuts", "", "#", BinSettings(15),"h_Cuts");

    t.CreateBranches(HistFac.makeTTree("tree"));

    kinfitter.SetZVertexSigma(0.0);
    LOG(INFO) << "Fit Z vertex enabled with sigma=0";

}

void EtapEPT::ProcessEvent(const TEvent& event, manager_t&)
{

    const TEventData& data = event.Reconstructed();

    h_Cuts->Fill("Seen",1.0);

    // start now with some cuts
    t.CBSumE = data.Trigger.CBEnergySum;
    t.CBAvgTime = data.Trigger.CBTiming;

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
    h_Cuts->Fill("1 in TAPS",1.0);

    std::sort(candidates.begin(), candidates.end(),
              [] (const TCandidatePtr& a, const TCandidatePtr& b) {
        return a->CaloEnergy > b->CaloEnergy;
    });

    // fill particles
    Particles_t particles;

    if(!findParticles(candidates, 2, particles, t, h_Cuts))
        return;

    // sum up the PID energy
    // (might be different to matched CB/PID Veto energy)
    t.PIDSumE = 0;
    for(const TCluster& cl : data.Clusters) {
        if(cl.DetectorType == Detector_t::Type_t::PID) {
            t.PIDSumE += cl.Energy;
        }
    }


    if(mc_smear && data.ID.isSet(TID::Flags_t::MC)) {
        particles.Proton = mc_smear->Smear(particles.Proton);
        for(auto& p : particles.Photons)
            p = mc_smear->Smear(p);
    }


    for(const TTaggerHit& taggerhit : data.TaggerHits) {
        promptrandom.SetTaggerHit(taggerhit.Time-t.CBAvgTime);
        promptrandom_wide.SetTaggerHit(taggerhit.Time);
        if(promptrandom.State() == PromptRandom::Case::Outside &&
           promptrandom_wide.State() == PromptRandom::Case::Outside)
            continue;

        t.TaggW = promptrandom.FillWeight();
        t.TaggW_wide = promptrandom_wide.FillWeight();
        t.TaggE = taggerhit.PhotonEnergy;
        t.TaggT = taggerhit.Time;
        t.TaggCh = taggerhit.Channel;


        for(unsigned ch=0;ch<EPT->GetNChannels();ch++) {
            if(noTaggChPerm && ch != t.TaggCh())
                continue;
            TTaggerHit taggerhit_(ch, EPT->GetPhotonEnergy(ch), taggerhit.Time);
            t.TaggCh_ = taggerhit_.Channel;
            t.TaggE_ = taggerhit_.PhotonEnergy;

            if(doKinfit(taggerhit_, kinfitter, particles, t, h_Cuts)) {
                t.IM_2g = particles.FittedPhotonSum.M();
                t.Tree->Fill();
            }
        }
    }
}

bool EtapEPT::findParticles(const TCandidatePtrList& candidates,
                               unsigned nPhotons,
                               EtapEPT::Particles_t& particles,
                               EtapEPT::Tree_t& t,
                               TH1D* h_CommonCuts)
{
    t.ProtonE = std_ext::NaN;
    t.ProtonTheta = std_ext::NaN;
    t.ProtonVetoE = std_ext::NaN;
    t.ProtonShortE = std_ext::NaN;
    t.DiscardedEk = std_ext::NaN;
    t.PhotonsEk = std_ext::NaN;
    t.nPhotonsCB   = 0;
    t.nPhotonsTAPS = 0;
    t.CBSumVetoE = std_ext::NaN;
    t.PhotonSum  = std_ext::NaN;
    t.PhotonThetas().resize(0);
    t.ProtonCopl = std_ext::NaN;

    h_CommonCuts->Fill("Seen Particles", 1.0);

    t.nCandidates = candidates.size();
    if(t.nCandidates<nPhotons+1)
        return false;
    h_CommonCuts->Fill("nCands ok", 1.0);
    h_CommonCuts->Fill("nCands exact", t.nCandidates==(nPhotons+1));

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
    photon_sum = {{0,0,0},0};
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
    const interval<double> ProtonCopl_cut(-45, 45);
    if(!ProtonCopl_cut.Contains(t.ProtonCopl))
        return false;
    h_CommonCuts->Fill("ProtonCopl ok", 1.0);

    return true;
}

bool EtapEPT::doKinfit(const TTaggerHit& taggerhit,
                          utils::KinFitter& kinfitter,
                          EtapEPT::Particles_t& particles,
                          EtapEPT::Tree_t& t,
                          TH1D* h_CommonCuts)
{
    h_CommonCuts->Fill("Seen KinFit", 1.0);

    // missing mass
    const LorentzVec beam_target = taggerhit.GetPhotonBeam() + LorentzVec({0, 0, 0}, ParticleTypeDatabase::Proton.Mass());
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
    h_CommonCuts->Fill("MM ok", 1.0);

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
    particles.FittedPhotonSum = {{0,0,0},0};

    for(const auto& photon : particles.FittedPhotons)
        particles.FittedPhotonSum += *photon;

    h_CommonCuts->Fill("KinFit ok", 1.0);

    return true;
}

void EtapEPT::ShowResult()
{
    canvas("Overview")
            << h_Cuts
            << TTree_drawable(t.Tree, "IM_2g >> h1(200,800,1050)","")
            << TTree_drawable(t.Tree, "IM_2g >> h2(200,800,1050)","TaggCh==TaggCh_")
            << TTree_drawable(t.Tree, "IM_2g >> h3(200,800,1050)","TaggW*(KinFitProb>0.01)")
            << TTree_drawable(t.Tree, "IM_2g >> h4(200,800,1050)","TaggW*(KinFitProb>0.01 && TaggCh == TaggCh_)")
            << TTree_drawable(t.Tree, "IM_2g >> h5(200,800,1050)","TaggW_wide*(KinFitProb>0.01 && TaggCh == TaggCh_)")
            << endc;

}

AUTO_REGISTER_PHYSICS(EtapEPT)
