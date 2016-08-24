#include "PullsWriter.h"

#include "analysis/plot/HistogramFactories.h"
#include "analysis/utils/Fitter.h"

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::utils;
using namespace std;

PullsWriter::PullTree_t& PullsWriter::getPullTree(const Fitter::FitParticle& particle)
{
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

PullsWriter::PullsWriter(HistogramFactory& histfac)
{
    photons_cb.CreateBranches(  histfac.makeTTree("pulls_photon_cb"));
    photons_taps.CreateBranches(histfac.makeTTree("pulls_photon_taps"));
    proton_cb.CreateBranches(   histfac.makeTTree("pulls_proton_cb"));
    proton_taps.CreateBranches( histfac.makeTTree("pulls_proton_taps"));

}

PullsWriter::~PullsWriter()
{}

void PullsWriter::Fill(const std::vector<Fitter::FitParticle>& fitParticles,
                       double tagger_weight, double fitprob)
{
    Fill(fitParticles, {}, tagger_weight, fitprob);
}

void PullsWriter::Fill(const std::vector<Fitter::FitParticle>& fitParticles,
                       const smear_sigmas_t& smear_sigmas,
                       double tagger_weight, double fitprob)
{

    const Fitter::FitParticle& proton = fitParticles.front();
    if(proton.Particle->Type() != ParticleTypeDatabase::Proton)
        throw std::runtime_error("First particle given is not a proton");

    for(const Fitter::FitParticle& p : fitParticles) {

        auto& tree = getPullTree(p);

        tree.TaggW = tagger_weight;
        tree.FitProb = fitprob;

//        tree.E     = p.Ek.Value_before;
//        tree.Theta = p.Theta.Value_before;
//        tree.Phi   = p.Phi.Value_before;
//        tree.ProtonTheta =  proton.Theta.Value_before;
//        tree.ProtonE = proton.Ek.Value_before;
        tree.ProtonTime = proton.Particle->Candidate->Time;

        const auto it_smear = smear_sigmas.find(p.Particle);
        if(it_smear == smear_sigmas.end()) {
//            tree.SigmaE = p.Ek.Sigma_before;
//            tree.SigmaTheta = p.Theta.Sigma_before;
//            tree.SigmaPhi = p.Phi.Sigma_before;
        }
        else {
            const Uncertainties_t& u = it_smear->second;
            tree.SigmaE = u.sigmaE;
            tree.SigmaTheta = u.sigmaTheta;
            tree.SigmaPhi = u.sigmaPhi;
        }

//        tree.PullE     = p.Ek.Pull;
//        tree.PullTheta = p.Theta.Pull;
//        tree.PullPhi   = p.Phi.Pull;

        tree.Multiplicity = fitParticles.size();

        tree.Tree->Fill();

    }
}
