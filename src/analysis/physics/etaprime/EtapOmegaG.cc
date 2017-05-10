#include "EtapOmegaG.h"

#include "utils/ParticleTools.h"
#include "utils/Matcher.h"
#include "utils/Combinatorics.h"
#include "analysis/utils/uncertainties/FitterSergey.h"
#include "analysis/utils/uncertainties/Interpolated.h"

#include "expconfig/ExpConfig.h"

#include "base/Logger.h"
#include "base/std_ext/math.h"
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
    APLCON::Fit_Settings_t settings;
    settings.MaxIterations = max_iterations;
    //    settings.ConstraintAccuracy = 1.0e-3;
    //    settings.Chi2Accuracy = 1.0e-2;
    //    settings.DebugLevel = 5;
    return settings;
}

EtapOmegaG::EtapOmegaG(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    promptrandom(ExpConfig::Setup::Get()), // use setup for promptrandom windows
    fitmodel_data(// use Interpolated, based on Sergey's model
                  utils::UncertaintyModels::Interpolated::makeAndLoad(
                      utils::UncertaintyModels::Interpolated::Type_t::Data,
                      // use Sergey as starting point
                      make_shared<utils::UncertaintyModels::FitterSergey>()
                      )),
    fitmodel_mc(// use Interpolated, based on Sergey's model
                utils::UncertaintyModels::Interpolated::makeAndLoad(
                    utils::UncertaintyModels::Interpolated::Type_t::MC,
                    // use Sergey as starting point
                    make_shared<utils::UncertaintyModels::FitterSergey>()
                    )),
    fitparams(true, // flag to enable z vertex
              3.0 // Z_vertex_sigma, =0 means unmeasured
              ),
    Sig(HistogramFactory("Sig",HistFac,"Sig"), fitparams),
    Ref(HistogramFactory("Ref",HistFac,"Ref"), fitparams)
{
    if(fitparams.Fit_Z_vertex) {
        LOG(INFO) << "Fit Z vertex enabled with sigma=" << fitparams.Z_vertex_sigma;
    }

    h_Cuts = HistFac.makeTH1D("Cuts", "", "#", BinSettings(15),"h_Cuts");
    h_DiscardedPhotons = HistFac.makeTH2D("DiscardedPhotons",
                                          {"#theta / #circ", BinSettings(50,0,8)},
                                          {"E_{kin} / MeV",  BinSettings(50,0,500)},
                                          "h_DiscardedPhotons");

    h_LostPhotons_sig = HistFac.makeTH1D("Sig: LostPhotons", "#theta", "#", BinSettings(200,0,180),"h_LostPhotons_sig");
    h_LostPhotons_ref = HistFac.makeTH1D("Ref: LostPhotons", "#theta", "#", BinSettings(200,0,180),"h_LostPhotons_ref");

    t.CreateBranches(Sig.treeCommon);
    t.CreateBranches(Ref.treeCommon);
    t.Tree = nullptr; // prevent accidental misuse...

    // setup does never change, so set it once and for all
    if(std_ext::contains(ExpConfig::Setup::Get().GetName(), "2014_07"))
        t.BeamTime = 1;
    else if(std_ext::contains(ExpConfig::Setup::Get().GetName(), "2014_10"))
        t.BeamTime = 2;
    else if(std_ext::contains(ExpConfig::Setup::Get().GetName(), "2014_12"))
        t.BeamTime = 3;
}

void EtapOmegaG::ProcessEvent(const TEvent& event, manager_t&)
{
    if(!triggersimu.ProcessEvent(event))
        h_Cuts->Fill("Triggersimu failed", 1.0);

    // we start with some general candidate handling,
    // later we split into ref/sig analysis according to
    // number of photons

    const bool have_MCTrue = !event.MCTrue().ID.IsInvalid();

    const TEventData& data = event.Reconstructed();

    const bool is_MC = data.ID.isSet(TID::Flags_t::MC);

    h_Cuts->Fill("Seen",1.0);
    h_Cuts->Fill(ExpConfig::Setup::Get().GetName().c_str(), 1.0);

    auto& particletree = event.MCTrue().ParticleTree;

    if(is_MC) {
        // until here, no physics cuts were done (THIS IS IMPORTANT)
        // so we can fill this into our mcWeightingEtaPrime instances
        Sig.mcWeightingEtaPrime.SetParticleTree(particletree);
        Ref.mcWeightingEtaPrime.SetParticleTree(particletree);
    }

    // Count EtaPrimes in MC sample
    h_Cuts->Fill("MCTrue #eta'", 0); // ensure the bin is there...
    if(particletree) {
        // note: this might also match to g p -> eta' eta' p,
        // but this is kinematically forbidden
        if(utils::ParticleTools::FindParticle(ParticleTypeDatabase::EtaPrime, particletree, 1)) {
            h_Cuts->Fill("MCTrue #eta'", 1);
        }
    }

    // do some MCTrue identification (if available)
    t.MCTrue = 0; // indicate data by default
    t.MCTrueMissed = "";
    t.TrueZVertex = event.MCTrue().Target.Vertex.z; // NaN in case of data

    if(particletree) {
        // 1=Signal, 2=Reference, 9=MissedBkg, >=10 found in ptreeBackgrounds
        if(particletree->IsEqual(ptreeSignal, utils::ParticleTools::MatchByParticleName)) {
            t.MCTrue = 1;
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
                t.MCTrueMissed = utils::ParticleTools::GetDecayString(particletree);
                t.MCTrue = 9;
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
        auto h_cut = t.MCTrue == 1 ? Sig.h_Cuts : Ref.h_Cuts;
        auto h_lost = t.MCTrue == 1 ? h_LostPhotons_sig : h_LostPhotons_ref;
        h_cut->Fill("MCTrue seen", 1.0);
        bool photons_accepted = true;
        auto mctrue_particles = utils::ParticleTypeList::Make(event.MCTrue().ParticleTree);
        for(const TParticlePtr& p : mctrue_particles.Get(ParticleTypeDatabase::Photon)) {
            if(geometry.DetectorFromAngles(*p) == Detector_t::Any_t::None) {
                h_lost->Fill(std_ext::radian_to_degree(p->Theta()));
                photons_accepted = false;
            }
        }
        if(photons_accepted) {
            h_cut->Fill("MCTrue Photon ok", 1.0);
        }
        auto proton = mctrue_particles.Get(ParticleTypeDatabase::Photon).front();
        if(geometry.DetectorFromAngles(*proton) != Detector_t::Any_t::None)
            h_cut->Fill("MCTrue Proton ok", 1.0);
    }

    // start now with some cuts

    if(!triggersimu.HasTriggered())
        return;
    h_Cuts->Fill("Triggered",1.0);

    t.CBSumE = triggersimu.GetCBEnergySum();
    t.CBAvgTime = triggersimu.GetRefTiming();

    if(data.Candidates.size()<3)
        return;
    h_Cuts->Fill("nCands>=3", 1.0);

    // etaprime physics always has something in TAPS
    // (the proton, by the way)
    // but remember: Backgrounds such as pi0pi0 have their proton in CB most likely
    // so it shouldn't be assumed that one of the TAPS clusters is the proton
    // (at least in the AntiPi0Pi0/AntiPi0Eta fits)
    bool haveTAPS = false;
    for(const auto& cand : data.Candidates.get_iter()) {
        if(cand->Detector & Detector_t::Type_t::TAPS) {
            haveTAPS = true;
        }
    }
    if(!haveTAPS)
        return;
    h_Cuts->Fill("1 in TAPS", 1.0);

    // sum up the PID energy
    // (might be different to matched CB/PID Veto energy)
    t.PIDSumE = 0;
    for(const TCluster& cl : data.Clusters) {
        if(cl.DetectorType == Detector_t::Type_t::PID) {
            t.PIDSumE += cl.Energy;
        }
    }

    // this ensures the TParticlePtr (shared_ptr) are only made once
    // but do not allow photons with polar angle <10degree
    utils::ProtonPhotonCombs proton_photons(data.Candidates, [this] (particle_t& p) {
        auto it = p.Photons.begin();
        while(it != p.Photons.end()) {
            auto& photon = *it;
            if(std_ext::radian_to_degree(photon->Theta())<7) {
                h_DiscardedPhotons->Fill(std_ext::radian_to_degree(photon->Theta()), photon->Ek());
                it = p.Photons.erase(it);
            }
            else
                ++it;
        }
    });

    // some extra info to pass to Process methods
    params_t p;
    p.MCTrue = t.MCTrue;
    p.ParticleTree = particletree;

    // set uncertainty model (maybe a bit ugly implemented here)
    Sig.kinfitter.SetUncertaintyModel(is_MC ? fitmodel_mc : fitmodel_data);
    Sig.treefitter_Pi0Pi0.SetUncertaintyModel(is_MC ? fitmodel_mc : fitmodel_data);
    Sig.treefitter_Pi0Eta.SetUncertaintyModel(is_MC ? fitmodel_mc : fitmodel_data);
    Sig.Pi0.treefitter.SetUncertaintyModel(is_MC ? fitmodel_mc : fitmodel_data);
    Sig.OmegaPi0.treefitter.SetUncertaintyModel(is_MC ? fitmodel_mc : fitmodel_data);
    Ref.kinfitter.SetUncertaintyModel(is_MC ? fitmodel_mc : fitmodel_data);


    for(const TTaggerHit& taggerhit : data.TaggerHits) {
        promptrandom.SetTaggerTime(triggersimu.GetCorrectedTaggerTime(taggerhit));
        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        t.TaggW  = promptrandom.FillWeight();
        t.TaggE  = taggerhit.PhotonEnergy;
        t.TaggT  = taggerhit.Time;
        t.TaggCh = taggerhit.Channel;
        t.TaggTcorr = triggersimu.GetCorrectedTaggerTime(taggerhit);

        p.TaggerHit = taggerhit;
        p.TaggW = t.TaggW;
        p.Particles = proton_photons(); // copy from pre-built combinations

        Sig.Process(p);
        Ref.Process(p);
    }

}

void EtapOmegaG::ProtonPhotonTree_t::Fill(const EtapOmegaG::params_t& params, const EtapOmegaG::particle_t& p, double fitted_proton_E)
{
    PhotonsEk = 0;
    nPhotonsCB = 0;
    nPhotonsTAPS = 0;
    CBSumVetoE = 0;
    PhotonThetas().clear();
    nTouchesHole = 0;
    for(const auto& photon : p.Photons) {
        const auto& cand = photon->Candidate;
        PhotonsEk += cand->CaloEnergy;
        if(cand->Detector & Detector_t::Type_t::CB) {
            nPhotonsCB++;
            CBSumVetoE += cand->VetoEnergy;
        }
        if(cand->Detector & Detector_t::Type_t::TAPS)
            nPhotonsTAPS++;
        PhotonThetas().emplace_back(std_ext::radian_to_degree(cand->Theta));
        nTouchesHole += cand->FindCaloCluster()->HasFlag(TCluster::Flags_t::TouchesHoleCentral);
    }
    assert(PhotonThetas().size() == p.Photons.size());

    DiscardedEk = p.DiscardedEk;
    PhotonSum = p.PhotonSum.M();
    MissingMass = p.MissingMass;
    ProtonCopl = std_ext::radian_to_degree(vec2::Phi_mpi_pi(p.Proton->Phi() - p.PhotonSum.Phi() - M_PI ));

    ProtonTime = p.Proton->Candidate->Time;
    ProtonE = p.Proton->Ek();
    ProtonTheta = std_ext::radian_to_degree(p.Proton->Theta());
    ProtonVetoE = p.Proton->Candidate->VetoEnergy;
    ProtonShortE = p.Proton->Candidate->FindCaloCluster()->ShortEnergy;
    auto true_proton = utils::ParticleTools::FindParticle(ParticleTypeDatabase::Proton, params.ParticleTree);
    if(true_proton)
        ProtonTrueAngle = std_ext::radian_to_degree(p.Proton->Angle(*true_proton));
    else
        ProtonTrueAngle = std_ext::NaN;

    FittedProtonE = fitted_proton_E;
}

EtapOmegaG::Sig_t::Sig_t(const HistogramFactory& HistFac, fitparams_t params) :
    h_Cuts(HistFac.makeTH1D("Cuts", "", "#", BinSettings(15),"h_Cuts")),
    h_MissedBkg(HistFac.makeTH1D("Missed Background", "", "", BinSettings(25),"h_MissedBkg")),
    treeCommon(HistFac.makeTTree("Common")),
    Pi0(params),
    OmegaPi0(params),
    mcWeightingEtaPrime(HistFac, utils::MCWeighting::EtaPrime),
    kinfitter(nullptr, params.Fit_Z_vertex,
              EtapOmegaG::MakeFitSettings(10)
              ),
    treefitter_Pi0Pi0(ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::TwoPi0_4g),
                      nullptr, params.Fit_Z_vertex, {},
                      MakeFitSettings(10)
                      ),
    treefitter_Pi0Eta(ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0Eta_4g),
                      nullptr, params.Fit_Z_vertex, {},
                      MakeFitSettings(10)
                      )
{
    const AxisSettings axis_IM("IM / MeV",{500,0,1000});
    h_IM_2g = HistFac.makeTH1D("IM 2#gamma",axis_IM,"h_IM_2g");
    h_IM_3g = HistFac.makeTH1D("IM 3#gamma",axis_IM,"h_IM_3g");
    h_IM_4g = HistFac.makeTH1D("IM 4#gamma",axis_IM,"h_IM_4g");

    t.CreateBranches(HistFac.makeTTree("Shared"));
    OmegaPi0.t.CreateBranches(HistFac.makeTTree("OmegaPi0"));
    Pi0.t.CreateBranches(HistFac.makeTTree("Pi0"));

    if(params.Fit_Z_vertex) {
        kinfitter.SetZVertexSigma(params.Z_vertex_sigma);
        treefitter_Pi0Pi0.SetZVertexSigma(params.Z_vertex_sigma);
        treefitter_Pi0Eta.SetZVertexSigma(params.Z_vertex_sigma);
    }

    {
        auto pi0s = treefitter_Pi0Pi0.GetTreeNodes(ParticleTypeDatabase::Pi0);
        treefitter_Pi0Pi0.SetIterationFilter([pi0s] () {
            auto lvsum1 = pi0s.front()->Get().LVSum;
            auto lvsum2 = pi0s.back()->Get().LVSum;

            const auto& pi0_cut = ParticleTypeDatabase::Pi0.GetWindow(80);

            return pi0_cut.Contains(lvsum1.M()) && pi0_cut.Contains(lvsum2.M());
        });
    }

    {
        auto pi0 = treefitter_Pi0Eta.GetTreeNode(ParticleTypeDatabase::Pi0);
        auto eta = treefitter_Pi0Eta.GetTreeNode(ParticleTypeDatabase::Eta);

        treefitter_Pi0Eta.SetIterationFilter([pi0,eta] () {
            const auto& pi0_lvsum = pi0->Get().LVSum;
            const auto& eta_lvsum = eta->Get().LVSum;

            const auto& pi0_cut = ParticleTypeDatabase::Pi0.GetWindow(80);
            const auto& eta_cut = ParticleTypeDatabase::Eta.GetWindow(120);

            return pi0_cut.Contains(pi0_lvsum.M()) && eta_cut.Contains(eta_lvsum.M());
        });
    }

}

void EtapOmegaG::Sig_t::Process(params_t params)
{
    params.Particles
            .Observe([this] (const std::string& s) { h_Cuts->Fill(s.c_str(), 1.0); }, "S ")
            .FilterMult(4, 70.0)
            .FilterIM({550, std_ext::inf})
            .FilterMM(params.TaggerHit, ParticleTypeDatabase::Proton.GetWindow(350).Round());

    if(params.Particles.empty())
        return;

    t.KinFitProb = std_ext::NaN;
    vector<double> IM_2g(6, std_ext::NaN);
    vector<double> IM_3g(4, std_ext::NaN);
    vector<double> IM_4g(1, std_ext::NaN);

    for(auto& p : params.Particles) {

        auto result = kinfitter.DoFit(params.TaggerHit.PhotonEnergy, p.Proton, p.Photons);

        if(result.Status != APLCON::Result_Status_t::Success)
            continue;

        if(!std_ext::copy_if_greater(t.KinFitProb, result.Probability))
            continue;

        t.KinFitProb = result.Probability;
        t.KinFitIterations = result.NIterations;
        t.KinFitZVertex = kinfitter.GetFittedZVertex();

        const auto& photons = kinfitter.GetFittedPhotons();
        utils::ParticleTools::FillIMCombinations(IM_2g.begin(), 2, photons);
        utils::ParticleTools::FillIMCombinations(IM_3g.begin(), 3, photons);
        utils::ParticleTools::FillIMCombinations(IM_4g.begin(), 4, photons);
    }

    if(!(t.KinFitProb > 0.005))
        return;

    for(auto& v : IM_2g)
        h_IM_2g->Fill(v, params.TaggW);
    for(auto& v : IM_3g)
        h_IM_3g->Fill(v, params.TaggW);
    for(auto& v : IM_4g)
        h_IM_4g->Fill(v, params.TaggW);

    h_Cuts->Fill("KinFit ok", 1.0);

    DoAntiPi0Eta(params);

    if(t.AntiPi0FitProb > 0.05)
        return;
    if(t.AntiEtaFitProb > 0.05)
        return;
    h_Cuts->Fill("Anti ok", 1.0);

    Pi0.Process(params);
    OmegaPi0.Process(params);

    if(!isfinite(Pi0.t.TreeFitProb) && !isfinite(OmegaPi0.t.TreeFitProb))
        return;

    h_Cuts->Fill("Sig ok", 1.0);

    if(isfinite(Pi0.t.TreeFitProb) && isfinite(OmegaPi0.t.TreeFitProb))
        h_Cuts->Fill("Both ok", 1.0);

    h_Cuts->Fill("Pi0 ok", isfinite(Pi0.t.TreeFitProb));
    h_Cuts->Fill("OmegaPi0 ok", isfinite(OmegaPi0.t.TreeFitProb));

    if(params.ParticleTree && params.MCTrue == 9) {
        const auto& decaystr = utils::ParticleTools::GetDecayString(params.ParticleTree);
        h_MissedBkg->Fill(decaystr.c_str(), 1.0);
    }

    // fill them all to keep them in sync
    treeCommon->Fill();
    t.Tree->Fill();
    Pi0.t.Tree->Fill();
    OmegaPi0.t.Tree->Fill();
    mcWeightingEtaPrime.Fill();
}

void EtapOmegaG::Sig_t::DoAntiPi0Eta(const params_t& params)
{
    t.AntiPi0FitProb = std_ext::NaN;
    t.AntiEtaFitProb = std_ext::NaN;

    for(const auto& p : params.Particles) {

        APLCON::Result_t r;

        treefitter_Pi0Pi0.PrepareFits(params.TaggerHit.PhotonEnergy,
                                      p.Proton, p.Photons);
        while(treefitter_Pi0Pi0.NextFit(r)) {
            if(r.Status != APLCON::Result_Status_t::Success)
                continue;
            if(!std_ext::copy_if_greater(t.AntiPi0FitProb, r.Probability))
                continue;
            // found fit with better prob
            t.AntiPi0FitIterations = r.NIterations;

            const auto& fitter = treefitter_Pi0Pi0;
            t.AntiPi0FitZVertex = fitter.GetFittedZVertex();
        }

        treefitter_Pi0Eta.PrepareFits(params.TaggerHit.PhotonEnergy,
                                      p.Proton, p.Photons);
        while(treefitter_Pi0Eta.NextFit(r)) {
            if(r.Status != APLCON::Result_Status_t::Success)
                continue;
            if(!std_ext::copy_if_greater(t.AntiEtaFitProb, r.Probability))
                continue;
            // found fit with better probability
            t.AntiEtaFitIterations = r.NIterations;

            const auto& fitter = treefitter_Pi0Eta;
            t.AntiEtaFitZVertex = fitter.GetFittedZVertex();
        }
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

    {
        treefitter.SetIterationFilter([this] () {
            const auto& pi0 = fitted_Pi0->Get().LVSum;
            double invchi2 = 1.0/std_ext::sqr(ParticleTypeDatabase::Pi0.Mass() - pi0.M());
            if(fitted_Omega) {
                const auto& omega = fitted_Omega->Get().LVSum;
                invchi2 += 1.0/std_ext::sqr(ParticleTypeDatabase::Omega.Mass() - omega.M());
            }
            return invchi2;
        },
        4);
    }
}

utils::TreeFitter EtapOmegaG::Sig_t::Fit_t::Make(const ParticleTypeDatabase::Type& subtree, fitparams_t params)
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
                EtapOmegaG::ptreeSignal,
                nullptr,
                params.Fit_Z_vertex,
                setupnodes,
                MakeFitSettings(15)
    };
    if(params.Fit_Z_vertex)
        treefitter.SetZVertexSigma(params.Z_vertex_sigma);
    return treefitter;
}

void fill_gNonPi0(
        EtapOmegaG::Sig_t::Fit_t::BaseTree_t& t,
        const TCandidatePtr& cand1, const TCandidatePtr& cand2)
{
    /// \todo on next round, change this to degree instead of radians
    t.gNonPi0_Theta().front() = cand1->Theta;
    t.gNonPi0_Theta().back()  = cand2->Theta;

    t.gNonPi0_CaloE().front() = cand1->CaloEnergy;
    t.gNonPi0_CaloE().back() = cand2->CaloEnergy;

    t.gNonPi0_VetoE().front() = cand1->VetoEnergy;
    t.gNonPi0_VetoE().back() = cand2->VetoEnergy;

    t.gNonPi0_TouchesHole().front() = cand1->FindCaloCluster()->HasFlag(TCluster::Flags_t::TouchesHoleCentral);
    t.gNonPi0_TouchesHole().back() = cand2->FindCaloCluster()->HasFlag(TCluster::Flags_t::TouchesHoleCentral);
}

void fill_PhotonCombs(EtapOmegaG::Sig_t::Fit_t::BaseTree_t& t, const TParticleList& photons)
{
    //  ggg combinatorics
    auto it_ggg = t.ggg().begin();
    for( auto comb = utils::makeCombination(photons,3); !comb.done(); ++comb ) {
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


EtapOmegaG::Sig_t::Pi0_t::Pi0_t(fitparams_t params) :
    Fit_t(Fit_t::Make(ParticleTypeDatabase::Pi0, params))
{

}

void EtapOmegaG::Sig_t::Pi0_t::Process(const params_t& params)
{
    t.TreeFitProb = std_ext::NaN;
    t.MCTrueMatch = 0;

    // for MCtrue identification
    TParticlePtr g1_Pi0_best;
    TParticlePtr g2_Pi0_best;
    TParticleList photons_best;

    for(const auto& p : params.Particles) {

        // do treefit
        treefitter.PrepareFits(params.TaggerHit.PhotonEnergy, p.Proton, p.Photons);

        APLCON::Result_t r;

        while(treefitter.NextFit(r)) {
            if(r.Status != APLCON::Result_Status_t::Success)
                continue;
            if(!std_ext::copy_if_greater(t.TreeFitProb, r.Probability))
                continue;
            // found fit with better prob
            t.TreeFitIterations = r.NIterations;
            t.TreeFitZVertex = treefitter.GetFittedZVertex();

            // for MCTrue matching
            g1_Pi0_best = fitted_g1_Pi0->Get().Leaf->Particle;
            g2_Pi0_best = fitted_g2_Pi0->Get().Leaf->Particle;
            photons_best = p.Photons;

            // IM fitted expected to be delta peaks since they were fitted...
            const LorentzVec& Pi0 = fitted_Pi0->Get().LVSum;
            t.IM_Pi0 = Pi0.M();

            t.IM_Pi0gg = fitted_EtaPrime->Get().LVSum.M();

            // there are two photon combinations possible for the omega
            // MC shows that it's the one with the higher IM_3g = IM_Pi0g
            auto leave1 = fitted_g_Omega->Get().Leaf;
            auto leave2 = fitted_g_EtaPrime->Get().Leaf;
            LorentzVec g1 = *leave1->AsFitted();
            LorentzVec g2 = *leave2->AsFitted();

            // invariant under swap
            t.IM_gg = (g1 + g2).M();
            const LorentzVec EtaPrime = g1 + g2 + Pi0;

            t.IM_Pi0g().front() = (Pi0 + g1).M();
            t.IM_Pi0g().back()  = (Pi0 + g2).M();
            if(t.IM_Pi0g().front() > t.IM_Pi0g().back()) {
                std::swap(t.IM_Pi0g().front(), t.IM_Pi0g().back());
                std::swap(leave1, leave2);
                std::swap(g1, g2);
            }

            // g1/leave1 is now the EtaPrime, g2/leave2 is now the Omega bachelor photon

            t.Bachelor_E().front() = Boost(g1, -EtaPrime.BoostVector()).E;
            t.Bachelor_E().back() =  Boost(g2, -EtaPrime.BoostVector()).E;

            fill_gNonPi0(t, leave1->Particle->Candidate, leave2->Particle->Candidate);
            fill_PhotonCombs(t, p.Photons);
            t.Fill(params, p, treefitter.GetFittedProton()->Ek());
        }

    }

    // there was at least one successful fit
    if(isfinite(t.TreeFitProb)) {

        // check MC matching
        if(params.MCTrue == 1) {
            auto& ptree_sig = params.ParticleTree;
            auto true_photons = utils::ParticleTools::FindParticles(ParticleTypeDatabase::Photon, ptree_sig);
            assert(true_photons.size() == 4);
            auto match_bycandidate = [] (const TParticlePtr& mctrue, const TParticlePtr& recon) {
                return mctrue->Angle(*recon->Candidate); // TCandidate converts into vec3
            };
            auto matched = utils::match1to1(true_photons, photons_best,
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


EtapOmegaG::Sig_t::OmegaPi0_t::OmegaPi0_t(fitparams_t params) :
    Fit_t(Fit_t::Make(ParticleTypeDatabase::Omega, params))
{

}

void EtapOmegaG::Sig_t::OmegaPi0_t::Process(const params_t& params)
{
    t.TreeFitProb = std_ext::NaN;
    t.MCTrueMatch = 0;

    // for MCTrue identification
    TParticlePtr g_Omega_best;
    TParticlePtr g_EtaPrime_best;
    TParticleList photons_best;

    // the EtaPrime bachelor photon is most important to us...
    TParticlePtr g_EtaPrime_fitted;
    LorentzVec EtaPrime_fitted;

    for(const auto& p : params.Particles) {

        // do treefit
        treefitter.PrepareFits(params.TaggerHit.PhotonEnergy, p.Proton, p.Photons);

        APLCON::Result_t r;

        while(treefitter.NextFit(r)) {
            if(r.Status != APLCON::Result_Status_t::Success)
                continue;
            if(!std_ext::copy_if_greater(t.TreeFitProb, r.Probability))
                continue;
            // found fit with better prob
            t.TreeFitIterations = r.NIterations;
            t.TreeFitZVertex = treefitter.GetFittedZVertex();

            // IM fitted expected to be delta peaks since they were fitted...
            EtaPrime_fitted = fitted_EtaPrime->Get().LVSum;
            t.IM_Pi0gg = EtaPrime_fitted.M();
            t.IM_Pi0g = fitted_Omega->Get().LVSum.M();
            t.IM_Pi0 = fitted_Pi0->Get().LVSum.M();

            // remember for matching
            g_EtaPrime_best = fitted_g_EtaPrime->Get().Leaf->Particle; // unfitted for matching
            g_Omega_best    = fitted_g_Omega->Get().Leaf->Particle; // unfitted for matching
            photons_best    = p.Photons;

            // have a look at the EtaPrime bachelor photon
            // the element NOT in the combination is the Bachelor photon
            g_EtaPrime_fitted = fitted_g_EtaPrime->Get().Leaf->AsFitted();

            t.IM_gg = ( *fitted_g_EtaPrime->Get().Leaf->AsFitted()
                        + *fitted_g_Omega->Get().Leaf->AsFitted()).M();

            fill_gNonPi0(t,
                         fitted_g_EtaPrime->Get().Leaf->Particle->Candidate,
                         fitted_g_Omega->Get().Leaf->Particle->Candidate);
            fill_PhotonCombs(t, p.Photons);
            t.Fill(params, p, treefitter.GetFittedProton()->Ek());
        }

    }

    // there was at least one successful fit
    if(isfinite(t.TreeFitProb)) {

        t.Bachelor_E = Boost(*g_EtaPrime_fitted, -EtaPrime_fitted.BoostVector()).E;

        // check MC matching
        if(params.MCTrue == 1) {
            auto& ptree_sig = params.ParticleTree;
            auto true_photons = utils::ParticleTools::FindParticles(ParticleTypeDatabase::Photon, ptree_sig);
            assert(true_photons.size() == 4);
            auto match_bycandidate = [] (const TParticlePtr& mctrue, const TParticlePtr& recon) {
                return mctrue->Angle(*recon->Candidate); // TCandidate converts into vec3
            };
            auto matched = utils::match1to1(true_photons, photons_best,
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
                auto g_EtaPrime = select_daughter(etap, ParticleTypeDatabase::Photon);
                auto omega = select_daughter(etap, ParticleTypeDatabase::Omega);
                auto g_Omega = select_daughter(omega, ParticleTypeDatabase::Photon);

                auto g_EtaPrime_matched = utils::FindMatched(matched, g_EtaPrime->Get());
                auto g_Omega_matched = utils::FindMatched(matched, g_Omega->Get());
                if(g_EtaPrime_matched == g_EtaPrime_best)
                    t.MCTrueMatch += 1;
                if(g_Omega_matched == g_Omega_best)
                    t.MCTrueMatch += 2;
            }
        }
    }
}

EtapOmegaG::Ref_t::Ref_t(const HistogramFactory& HistFac, EtapOmegaG::fitparams_t params) :
    h_Cuts(HistFac.makeTH1D("Cuts", "", "#", BinSettings(15),"h_Cuts")),
    h_MissedBkg(HistFac.makeTH1D("Missed Background", "", "", BinSettings(25),"h_MissedBkg")),
    treeCommon(HistFac.makeTTree("Common")),
    mcWeightingEtaPrime(HistFac, utils::MCWeighting::EtaPrime),
    kinfitter(nullptr, params.Fit_Z_vertex,
              EtapOmegaG::MakeFitSettings(15)
              )
{
    t.CreateBranches(HistFac.makeTTree("Ref"));
    if(params.Fit_Z_vertex)
        kinfitter.SetZVertexSigma(params.Z_vertex_sigma);
}

void EtapOmegaG::Ref_t::Process(params_t params)
{
    params.Particles
            .Observe([this] (const std::string& s) { h_Cuts->Fill(s.c_str(), 1.0); }, "R ")
            .FilterMult(2, 70.0)
            .FilterIM({600, std_ext::inf})
            .FilterMM(params.TaggerHit, ParticleTypeDatabase::Proton.GetWindow(350).Round());

    if(params.Particles.empty())
        return;

    t.KinFitProb = std_ext::NaN;
    for(const auto& p : params.Particles) {

        auto result = kinfitter.DoFit(params.TaggerHit.PhotonEnergy, p.Proton, p.Photons);

        if(result.Status != APLCON::Result_Status_t::Success)
            continue;

        if(!std_ext::copy_if_greater(t.KinFitProb, result.Probability))
            continue;

        t.KinFitProb = result.Probability;
        t.KinFitIterations = result.NIterations;
        t.KinFitZVertex = kinfitter.GetFittedZVertex();

        t.Fill(params, p, kinfitter.GetFittedProton()->Ek());

        const auto& fittedPhotons = kinfitter.GetFittedPhotons();
        t.IM_2g = (*fittedPhotons.front() + *fittedPhotons.back()).M();
        t.IM_2g_raw = (*p.Photons.front() + *p.Photons.back()).M();
    }

    if(t.KinFitProb>0.005) {

        if(params.ParticleTree && params.MCTrue == 9) {
            const auto& decaystr = utils::ParticleTools::GetDecayString(params.ParticleTree);
            h_MissedBkg->Fill(decaystr.c_str(), 1.0);
        }

        h_Cuts->Fill("Fill", 1.0);
        treeCommon->Fill();
        t.Tree->Fill();
        mcWeightingEtaPrime.Fill();
    }

}

void EtapOmegaG::ShowResult()
{
    canvas(GetName()+": Overview")
            << h_Cuts << drawoption("colz") << h_DiscardedPhotons << endr
            << Sig.h_Cuts << Sig.h_MissedBkg << h_LostPhotons_sig
            << endr
            << Ref.h_Cuts << Ref.h_MissedBkg << h_LostPhotons_ref
            << endc;


    const string mcWeightRef = Ref.mcWeightingEtaPrime.FriendTTree(Ref.t.Tree) ? "MCWeight" : "";

    canvas(GetName()+": Reference")
            << TTree_drawable(Ref.t.Tree, "IM_2g >> (100,800,1050)", mcWeightRef)
            << TTree_drawable(Ref.t.Tree, "IM_2g_raw >> (100,800,1050)", mcWeightRef)
            << endc;

    Sig.Pi0.t.Tree->AddFriend(Sig.t.Tree);
    Sig.OmegaPi0.t.Tree->AddFriend(Sig.t.Tree);

    const string mcWeightSig = Sig.mcWeightingEtaPrime.FriendTTree(Sig.OmegaPi0.t.Tree) &&
                               Sig.mcWeightingEtaPrime.FriendTTree(Sig.Pi0.t.Tree)
                               ? "*MCWeight" : "";

    canvas(GetName()+": Signal")
            << TTree_drawable(Sig.OmegaPi0.t.Tree, "Bachelor_E >> (100,50,250)","(TreeFitProb>0.01)"+mcWeightSig)
            << TTree_drawable(Sig.Pi0.t.Tree, "Bachelor_E[0] >> (100,50,250)","(TreeFitProb>0.01)"+mcWeightSig)
            << endr
            << TTree_drawable(Sig.OmegaPi0.t.Tree, "IM_Pi0gg >> (100,800,1050)","(TreeFitProb>0.01)"+mcWeightSig)
            << TTree_drawable(Sig.Pi0.t.Tree, "IM_Pi0gg >> (100,800,1050)","(TreeFitProb>0.01)"+mcWeightSig)
            << endr
            << TTree_drawable(Sig.OmegaPi0.t.Tree, "MCTrueMatch")
            << TTree_drawable(Sig.Pi0.t.Tree, "MCTrueMatch")
            << endc;
}

void EtapOmegaG::Finish()
{
    Sig.mcWeightingEtaPrime.Finish();
    Ref.mcWeightingEtaPrime.Finish();
}


const std::vector<EtapOmegaG::Background_t> EtapOmegaG::ptreeBackgrounds = {
    {"1Pi0", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0_2g)},
    {"2Pi0", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::TwoPi0_4g)},
    {"Pi0Eta_4g", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0Eta_4g)},
    {"3Pi0", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::ThreePi0_6g)},
    {"OmegaPi0g", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Omega_gPi0_3g)},
    {"OmegaPi0PiPPiM", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Omega_Pi0PiPPiM_2g)},
    {"EtaP2Pi0Eta", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::EtaPrime_2Pi0Eta_6g)},
    {"2Pi0Dalitz", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::TwoPi0_2ggEpEm)},
    {"3Pi0Dalitz", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::ThreePi0_4ggEpEm)},
    {"1Eta", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Eta_2g)},
    {"Eta3Pi0", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Eta_3Pi0_6g)},
    {"Pi0Eta_8g", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0Eta_Pi03Pi0_8g)},
    {"EtaP2gPiPPiM", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::EtaPrime_EtaPiPPiM_2gPiPPiM)},
    {"Pi0Eta_eeg2g", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0Eta_gEpEm2g)},
    {"Pi0PiPPiM", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0PiPi_2gPiPi)},
    {"RhoPiPPiM", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Rho_PiPi)},
    {"2Pi0PiPPiM", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::TwoPi0PiPi_4gPiPi)},
};

AUTO_REGISTER_PHYSICS(EtapOmegaG)
