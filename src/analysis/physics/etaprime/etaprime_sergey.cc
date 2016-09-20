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
    Debug(opts->Get<bool>("Debug", false)),
    params(// use FitterSergey as default
           make_shared<utils::UncertaintyModels::FitterSergey>(),
           true, // flag to enable z vertex
           3.0 // Z_vertex_sigma, =0 means unmeasured
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

    treeSergey.CreateBranches(HistFac.makeTTree("treeSergey"));
    treeAnt.CreateBranches(HistFac.makeTTree("treeAnt"));

    h_MissedBkg = HistFac.makeTH1D("Missed Background", "", "#", BinSettings(10),"h_MissedBkg");
}

void EtapSergey::fillTree(EtapSergey::Tree_t& t,
                          const std::vector<EtapSergey::result_t>& results,
                          unsigned MCTrue)
{
    t.MCTrue = MCTrue;
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
        t.PIDSumE    = r.PIDSumE;

        t.Tree->Fill();
    }
}

void EtapSergey::ProcessEvent(const TEvent& event, manager_t&)
{
    const TEventData& data = event.Reconstructed();

    // run Sergey's analysis
    auto results_sergey = fitter_sergey.Process(data);

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

    fillTree(treeSergey, results_sergey, MCTrue);

    vector<result_t> results_ant;

    // try to reproduce the result from Sergey
    for(const auto& r_sergey : results_sergey) {

        result_t r;
        r.TaggE = r_sergey.TaggE;
        r.TaggCh = r_sergey.TaggCh;
        r.TaggT = r_sergey.TaggT;

        r.KinFitProb = std_ext::NaN;

        for(const auto& cand_proton :  data.Candidates.get_iter()) {

            auto proton = make_shared<TParticle>(ParticleTypeDatabase::Proton, cand_proton);

            TParticleList photons;
            for(const auto& cand : data.Candidates.get_iter()) {
                if(cand == cand_proton)
                    continue;
                auto photon = make_shared<TParticle>(ParticleTypeDatabase::Photon, cand);
                photons.emplace_back(photon);
            }

            kinfitter.SetEgammaBeam(r.TaggE);
            kinfitter.SetProton(proton);
            kinfitter.SetPhotons(photons);

            auto result = kinfitter.DoFit();

            if(!std_ext::copy_if_greater(r.KinFitProb, result.Probability))
                continue;

            const auto fitted_photons = kinfitter.GetFittedPhotons();
            LorentzVec fitted_photon_sum;
            for(auto& p : fitted_photons)
                fitted_photon_sum += *p;
            r.IM_4g = fitted_photon_sum.M();
        }

        results_ant.emplace_back(move(r));
    }

    fillTree(treeAnt, results_ant, MCTrue);

}

void EtapSergey::ShowResult()
{
    canvas("Result")
            << TTree_drawable(treeSergey.Tree, "IM_4g >> (100,800,1050)","KinFitProb>0.01")
            << TTree_drawable(treeAnt.Tree,    "IM_4g >> (100,800,1050)","KinFitProb>0.01")
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
