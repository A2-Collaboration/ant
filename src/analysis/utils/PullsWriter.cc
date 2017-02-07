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
                       double tagger_weight, double fitprob, double fitted_z_vertex)
{

    const TParticlePtr& proton = fitParticles.front().Particle;
    if(proton->Type() != ParticleTypeDatabase::Proton)
        throw std::runtime_error("First particle given is not a proton");

    for(const Fitter::FitParticle& p : fitParticles) {
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

        tree.ProtonE = proton->Ek();
        tree.ProtonTheta = proton->Theta();
        tree.ProtonTime = proton->Candidate->Time;

        tree.Values = p.GetValues();
        tree.Sigmas = p.GetSigmas();
        tree.Pulls = p.GetPulls();

        tree.Tree->Fill();

    }
}
