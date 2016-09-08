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
           3.0 // Z_vertex_sigma, =0 means unmeasured
           ),
    kinfitter_sig("kinfitter_sig",4,
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
        kinfitter_sig.SetZVertexSigma(params.Z_vertex_sigma);
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

}

void EtapSergey::ProcessEvent(const TEvent& event, manager_t&)
{
    const TEventData& data = event.Reconstructed();

    auto results_sergey = fitter_sergey.Process(data);

    if(results_sergey.empty())
        return;

    cout << ">>>> Sergey:" << endl;
    for(auto& r : results_sergey)
        cout << r << endl;

    vector<double> photon_energies(results_sergey.size());
    std::transform(results_sergey.begin(), results_sergey.end(), photon_energies.begin(),
                   [] (const utils::FitterSergey::result_t& r) {
        return r.TaggE;
    });

    vector<utils::FitterSergey::result_t> results_ant;

    // try to reproduce the result from Sergey
    for(const auto& photon_energy : photon_energies) {

        utils::FitterSergey::result_t r;

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

            kinfitter_sig.SetEgammaBeam(photon_energy);
            kinfitter_sig.SetProton(proton);
            kinfitter_sig.SetPhotons(photons);

            auto result = kinfitter_sig.DoFit();

            if(!std_ext::copy_if_greater(r.KinFitProb, result.Probability))
                continue;
        }

        results_ant.emplace_back(move(r));
    }


    cout << ">>>> Ant:" << endl;
    for(auto& r : results_ant)
        cout << r << endl;
    cout << endl;
}

void EtapSergey::ShowResult()
{

}

AUTO_REGISTER_PHYSICS(EtapSergey)
