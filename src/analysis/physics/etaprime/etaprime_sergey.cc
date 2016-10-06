#include "etaprime_sergey.h"

#include "plot/root_draw.h"
#include "base/std_ext/misc.h"
#include "base/Logger.h"
#include "utils/particle_tools.h"

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace std;

utils::TreeFitter Make(const EtapSergey::params_t& params)
{
    auto setupnodes = [] (const ParticleTypeTree& t) {
        utils::TreeFitter::nodesetup_t nodesetup;
        // always exlude the EtaPrime
        if(t->Get() == ParticleTypeDatabase::EtaPrime)
            nodesetup.Excluded = true;
        // exclude the Omega as well
        if(t->Get() == ParticleTypeDatabase::Omega)
            nodesetup.Excluded = true;
        return nodesetup;
    };

    utils::TreeFitter treefitter{
        "sig_treefitter_Pi0",
        ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::EtaPrime_gOmega_ggPi0_4g),
        params.Fit_uncertainty_model,
        params.Fit_Z_vertex,
        setupnodes
    };
    if(params.Fit_Z_vertex)
        treefitter.SetZVertexSigma(params.Z_vertex_sigma);
    return treefitter;
}



EtapSergey::EtapSergey(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    params(// use FitterSergey as default
           make_shared<utils::UncertaintyModels::FitterSergey>(),
           true, // flag to enable z vertex
           0.333*10.0 // Z_vertex_sigma, =0 means unmeasured
           ),
    kinfitter("kinfitter_sig",4,
                  params.Fit_uncertainty_model, params.Fit_Z_vertex
                  ),
    treefitter(Make(params)),
    treefitter_Pi0Pi0("treefit_Pi0Pi0",
                      ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::TwoPi0_4g),
                      params.Fit_uncertainty_model, params.Fit_Z_vertex, {}
                      ),
    treefitter_Pi0Eta("treefit_Pi0Eta",
                      ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0Eta_4g),
                      params.Fit_uncertainty_model, params.Fit_Z_vertex, {}
                      ),
    fitted_Pi0(treefitter.GetTreeNode(ParticleTypeDatabase::Pi0)),
    fitted_Omega(treefitter.GetTreeNode(ParticleTypeDatabase::Omega)),
    fitted_EtaPrime(treefitter.GetTreeNode(ParticleTypeDatabase::EtaPrime))
{
    if(params.Fit_Z_vertex) {
        kinfitter.SetZVertexSigma(params.Z_vertex_sigma);
        treefitter_Pi0Pi0.SetZVertexSigma(params.Z_vertex_sigma);
        treefitter_Pi0Eta.SetZVertexSigma(params.Z_vertex_sigma);
    }

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
            const interval<double> pi0_cut{75,210};
            return pi0_cut.Contains(pi0.M());
        });
    }

    {
        auto pi0s = treefitter_Pi0Pi0.GetTreeNodes(ParticleTypeDatabase::Pi0);
        treefitter_Pi0Pi0.SetIterationFilter([pi0s] () {
            auto lvsum1 = pi0s.front()->Get().LVSum;
            auto lvsum2 = pi0s.back()->Get().LVSum;

            const interval<double> pi0_cut{75,210};

            return pi0_cut.Contains(lvsum1.M()) && pi0_cut.Contains(lvsum2.M());
        });
    }

    {
        auto pi0 = treefitter_Pi0Eta.GetTreeNode(ParticleTypeDatabase::Pi0);
        auto eta = treefitter_Pi0Eta.GetTreeNode(ParticleTypeDatabase::Eta);

        treefitter_Pi0Eta.SetIterationFilter([pi0,eta] () {
            const auto& pi0_lvsum = pi0->Get().LVSum;
            const auto& eta_lvsum = eta->Get().LVSum;

            const interval<double> pi0_cut{75,210};
            const interval<double> eta_cut{370,675};

            return pi0_cut.Contains(pi0_lvsum.M()) && eta_cut.Contains(eta_lvsum.M());
        });
    }

    treeSergey.CreateBranches(HistFac.makeTTree("treeSergey"));
    treeAnt.CreateBranches(HistFac.makeTTree("treeAnt"));

    h_MissedBkg = HistFac.makeTH1D("Missed Background", "", "#", BinSettings(10),"h_MissedBkg");
}

void EtapSergey::fillTree(EtapSergey::Tree_t& t,
                          const std::vector<EtapSergey::result_t>& results,
                          unsigned MCTrue, double PIDSumE)
{
    t.MCTrue = MCTrue;
    t.PIDSumE = PIDSumE;

    for(const result_t& r : results) {
        t.TaggE  = r.TaggE;
        t.TaggT  = r.TaggT;
        t.TaggCh = r.TaggCh;

        t.KinFitProb     = r.KinFitProb;
        t.TreeFitProb    = r.TreeFitProb;
        t.AntiPi0FitProb = r.AntiPi0FitProb;
        t.AntiEtaFitProb = r.AntiEtaFitProb;

        t.IM_3g = r.IM_3g;
        t.IM_4g = r.IM_4g;

        t.gNonPi0_Theta = r.gNonPi0_Theta;
        t.gNonPi0_CaloE = r.gNonPi0_CaloE;

        t.CBVetoSumE = r.CBVetoSumE;

        t.KinFitProtonIdx = r.KinFitProtonIdx;
        t.TreeFitProtonIdx = r.TreeFitProtonIdx;
        t.AntiPi0FitProtonIdx = r.AntiPi0FitProtonIdx;
        t.AntiEtaFitProtonIdx = r.AntiEtaFitProtonIdx;

        t.Tree->Fill();
    }
}

void EtapSergey::ProcessEvent(const TEvent& event, manager_t&)
{
    const TEventData& data = event.Reconstructed();

    TCandidatePtrList cands;
    for(auto& cand : data.Candidates.get_iter()) {
        if(cand->Detector & Detector_t::Type_t::CB)
            cands.push_back(cand);
    }
    for(auto& cand : data.Candidates.get_iter()) {
        if(cand->Detector & Detector_t::Type_t::TAPS)
            cands.push_back(cand);
    }

    // run Sergey's analysis
    auto results_sergey = fitter_sergey.Process(data.TaggerHits, cands);

    if(results_sergey.empty())
        return;

    // do some MCTrue identification (if available)
    auto MCTrue = 0; // indicate data by default
    {
        auto& particletree = event.MCTrue().ParticleTree;
        if(particletree) {
            // 1=Signal, 2=Reference, 9=MissedBkg, >=10 found in ptreeBackgrounds
            if(particletree->IsEqual(ptreeSignal, utils::ParticleTools::MatchByParticleName)) {
                MCTrue = 1;
            }
            else {
                MCTrue = 10;
                bool found = false;
                for(const auto& ptreeBkg : ptreeBackgrounds) {
                    if(particletree->IsEqual(ptreeBkg.Tree, utils::ParticleTools::MatchByParticleName)) {
                        found = true;
                        break;
                    }
                    MCTrue++;
                }
                if(!found) {
                    MCTrue = 9;
                    const auto& decaystr = utils::ParticleTools::GetDecayString(particletree);
                    h_MissedBkg->Fill(decaystr.c_str(), 1.0);
                }
            }
        }
        else if(!event.MCTrue().ID.IsInvalid()) {
            // in rare cases, the particletree is not available, although we're running on MCTrue
            // mark this as other MC background
            MCTrue = 9;
        }
    }


    double PIDSumE = 0;
    for(const TCluster& cl : data.Clusters) {
        if(cl.DetectorType == Detector_t::Type_t::PID) {
            PIDSumE += cl.Energy;
        }
    }

    fillTree(treeSergey, results_sergey,
             MCTrue, PIDSumE);

    vector<result_t> results_ant;

    // build permuted protons
    struct particles_t {
        particles_t(const TParticlePtr& proton, unsigned protonIdx) :
            Proton(proton), ProtonIdx(protonIdx) {}
        TParticleList Photons;
        TParticlePtr  Proton;
        unsigned      ProtonIdx;
    };

    std::vector<particles_t> particles;
    {
        TParticleList all_photons;
        TParticleList all_protons;

        for(const auto& cand_proton :  cands) {
            all_protons.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::Proton, cand_proton));
            all_photons.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::Photon, cand_proton));
        }

        for(auto it_proton = all_protons.cbegin(); it_proton != all_protons.cend(); ++it_proton) {
            auto& proton = *it_proton;
            particles.emplace_back(proton, std::distance(all_protons.cbegin(), it_proton)+1);
            auto& photons = particles.back().Photons;
            for(const auto& photon : all_photons) {
                if(proton->Candidate == photon->Candidate)
                    continue;
                photons.emplace_back(photon);
            }
        }
    }


    // try to reproduce the result from Sergey
    for(const auto& r_sergey : results_sergey) {

        result_t r;
        r.TaggE = r_sergey.TaggE;
        r.TaggCh = r_sergey.TaggCh;
        r.TaggT = r_sergey.TaggT;

        // KinFit
        r.KinFitProb = std_ext::NaN;
        for(const auto& p : particles) {
            kinfitter.SetEgammaBeam(r.TaggE);
            kinfitter.SetProton(p.Proton);
            kinfitter.SetPhotons(p.Photons);

            auto result = kinfitter.DoFit();
            if(result.Status != APLCON::Result_Status_t::Success)
                continue;
            if(!std_ext::copy_if_greater(r.KinFitProb, result.Probability))
                continue;

            r.KinFitProtonIdx = p.ProtonIdx;
        }

        // AntiPi0Fit
        r.AntiPi0FitProb = std_ext::NaN;
        for(const auto& p : particles) {
            treefitter_Pi0Pi0.SetEgammaBeam(r.TaggE);
            treefitter_Pi0Pi0.SetProton(p.Proton);
            treefitter_Pi0Pi0.SetPhotons(p.Photons);

            APLCON::Result_t result;
            while(treefitter_Pi0Pi0.NextFit(result)) {
                if(result.Status != APLCON::Result_Status_t::Success)
                    continue;
                if(!std_ext::copy_if_greater(r.AntiPi0FitProb, result.Probability))
                    continue;
                r.AntiPi0FitProtonIdx = p.ProtonIdx;
            }
        }

        // AntiEtaFit
        r.AntiEtaFitProb = std_ext::NaN;
        for(const auto& p : particles) {
            treefitter_Pi0Eta.SetEgammaBeam(r.TaggE);
            treefitter_Pi0Eta.SetProton(p.Proton);
            treefitter_Pi0Eta.SetPhotons(p.Photons);

            APLCON::Result_t result;
            while(treefitter_Pi0Eta.NextFit(result)) {
                if(result.Status != APLCON::Result_Status_t::Success)
                    continue;
                if(!std_ext::copy_if_greater(r.AntiEtaFitProb, result.Probability))
                    continue;
                r.AntiEtaFitProtonIdx = p.ProtonIdx;
            }
        }

        // TreeFit
        r.TreeFitProb = std_ext::NaN;
        for(const auto& p : particles) {
            treefitter.SetEgammaBeam(r.TaggE);
            treefitter.SetProton(p.Proton);
            treefitter.SetPhotons(p.Photons);

            APLCON::Result_t result;
            while(treefitter.NextFit(result)) {
                if(result.Status != APLCON::Result_Status_t::Success)
                    continue;
                if(!std_ext::copy_if_greater(r.TreeFitProb, result.Probability))
                    continue;
                r.TreeFitProtonIdx = p.ProtonIdx;
                const auto fitted_photons = treefitter.GetFittedPhotons();
                LorentzVec fitted_photon_sum;
                for(auto& p : fitted_photons)
                    fitted_photon_sum += *p;
                r.IM_4g = fitted_photon_sum.M();
            }
        }

        results_ant.emplace_back(move(r));
    }

    fillTree(treeAnt, results_ant, MCTrue, PIDSumE);

}

void EtapSergey::ShowResult()
{
    treeAnt.Tree->AddFriend(treeSergey.Tree);

    canvas("Result")
            << drawoption("colz")
            << TTree_drawable(treeSergey.Tree, "IM_4g >> (100,700,1050)","KinFitProb>0.01")
            << TTree_drawable(treeAnt.Tree, "treeSergey.IM_4g:treeAnt.IM_4g >> (100,800,1050,100,800,1050)","treeSergey.KinFitProb>0.01 && treeAnt.KinFitProb>0.01")
            << TTree_drawable(treeAnt.Tree, "treeSergey.KinFitProb:treeAnt.KinFitProb >> (100,0,1,100,0,1)","")
            << TTree_drawable(treeAnt.Tree, "treeSergey.KinFitProtonIdx:treeAnt.KinFitProtonIdx >> (5,1,6,5,1,6)","")
            << TTree_drawable(treeAnt.Tree, "treeSergey.AntiPi0FitProb:treeAnt.AntiPi0FitProb >> (100,0,0.002,100,0,0.002)","")
            << TTree_drawable(treeAnt.Tree, "treeSergey.AntiPi0FitProtonIdx:treeAnt.AntiPi0FitProtonIdx >> (5,1,6,5,1,6)","")
            << TTree_drawable(treeAnt.Tree, "treeSergey.AntiEtaFitProb:treeAnt.AntiEtaFitProb >> (100,0,0.002,100,0,0.002)","")
            << TTree_drawable(treeAnt.Tree, "treeSergey.AntiEtaFitProtonIdx:treeAnt.AntiEtaFitProtonIdx >> (5,1,6,5,1,6)","")
            << TTree_drawable(treeAnt.Tree, "treeSergey.TreeFitProb:treeAnt.TreeFitProb >> (100,0,1,100,0,1)","")
            << TTree_drawable(treeAnt.Tree, "treeSergey.TreeFitProtonIdx:treeAnt.TreeFitProtonIdx >> (5,1,6,5,1,6)","")
            << endc;
}

const ParticleTypeTree EtapSergey::ptreeSignal = ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::EtaPrime_gOmega_ggPi0_4g);

const std::vector<EtapSergey::Background_t> EtapSergey::ptreeBackgrounds = {
    {"1Pi0", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0_2g)},
    {"2Pi0", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::TwoPi0_4g)},
    {"Pi0Eta", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0Eta_4g)},
    {"3Pi0", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::ThreePi0_6g)},
    {"Pi0PiPi", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0PiPi_2gPiPi)},
    {"OmegaPi0g", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Omega_gPi0_3g)},
    {"OmegaPi0PiPPiM", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Omega_Pi0PiPPiM_2g)},
    {"EtaP2Pi0Eta", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::EtaPrime_2Pi0Eta_6g)},
    {"2Pi0Dalitz", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::TwoPi0_2ggEpEm)},
    {"3Pi0Dalitz", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::ThreePi0_4ggEpEm)},
    {"1Eta", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Eta_2g)},
};

AUTO_REGISTER_PHYSICS(EtapSergey)
