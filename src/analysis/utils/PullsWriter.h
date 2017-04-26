#pragma once

#include "base/WrapTTree.h"
#include "analysis/utils/fitter/Fitter.h"
#include "analysis/plot/HistogramFactory.h"

namespace ant {
namespace analysis {

class HistogramFactory;

namespace utils {

class PullsWriter {

public:
    struct PullTree_t : WrapTTree {

        ADD_BRANCH_T(double, TaggW) // for prompt-random subtraction
        ADD_BRANCH_T(double, FitProb)
        ADD_BRANCH_T(double, FittedZVertex)
        ADD_BRANCH_T(unsigned, Multiplicity)

        ADD_BRANCH_T(double, E)
        ADD_BRANCH_T(double, Theta)
        ADD_BRANCH_T(double, Phi)
        ADD_BRANCH_T(double, ShowerDepth)

        // see utils::Fitter::FitParticle for the parametrization
        ADD_BRANCH_T(std::vector<double>, Values)
        ADD_BRANCH_T(std::vector<double>, Sigmas)
        ADD_BRANCH_T(std::vector<double>, Pulls)
    };

protected:
    PullTree_t photons_cb;
    PullTree_t photons_taps;

    PullTree_t proton_cb;
    PullTree_t proton_taps;

    PullTree_t& getPullTree(const utils::Fitter::FitParticle& particle) {
        const auto& det = particle.Particle->Candidate->Detector;

        if(particle.Particle->Type() == ParticleTypeDatabase::Photon)
        {

            if(det & Detector_t::Type_t::CB) {
                return photons_cb;
            } else if(det & Detector_t::Type_t::TAPS) {
                return photons_taps;
            } else
                throw std::runtime_error("Unexpected detector type in fitter");

        }
        else if(particle.Particle->Type() == ParticleTypeDatabase::Proton)
        {

            if(det & Detector_t::Type_t::CB) {
                return proton_cb;
            } else if(det & Detector_t::Type_t::TAPS) {
                return proton_taps;
            } else
                throw std::runtime_error("Unexpected detector type in fitter");

        }
        else
            throw std::runtime_error("Unexpected Particle type in fitter!");
    }

public:

    PullsWriter(const HistogramFactory& histfac) {
        photons_cb.CreateBranches(  histfac.makeTTree("pulls_photon_cb"));
        photons_taps.CreateBranches(histfac.makeTTree("pulls_photon_taps"));
        proton_cb.CreateBranches(   histfac.makeTTree("pulls_proton_cb"));
        proton_taps.CreateBranches( histfac.makeTTree("pulls_proton_taps"));
    }

    void Fill(const std::vector<Fitter::FitParticle>& fitParticles,
              double tagger_weight, double fitprob, double fitted_z_vertex)
    {
        const auto& proton = fitParticles.front().Particle;
        if(proton->Type() != ParticleTypeDatabase::Proton)
            throw std::runtime_error("First particle given is not a proton");

        for(const auto& p : fitParticles) {
            if(p.Particle->Candidate->FindCaloCluster()->HasFlag(TCluster::Flags_t::TouchesHoleCentral))
               continue;

            auto& tree = getPullTree(p);

            tree.TaggW = tagger_weight;
            tree.FitProb = fitprob;
            tree.FittedZVertex = fitted_z_vertex;
            tree.Multiplicity = fitParticles.size();

            tree.E = p.Particle->Ek();
            tree.Theta = p.Particle->Theta();
            tree.Phi = p.Particle->Phi();
            tree.ShowerDepth = p.GetShowerDepth();

            tree.Values = p.GetValues_before();
            tree.Sigmas = p.GetSigmas_before();
            tree.Pulls = p.GetPulls();

            tree.Tree->Fill();

        }
    }

};

}
}
}
