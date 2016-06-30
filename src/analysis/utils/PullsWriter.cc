#include "PullsWriter.h"

#include "analysis/plot/HistogramFactories.h"
#include "analysis/utils/Fitter.h"

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::utils;
using namespace std;

//template <typename val, typename data>
//data& Sel(const val& V, const val& A, const val& B, data& dA, data& dB) {

//    if(V == A) {
//        return dA;
//    } else if(V == B) {
//        return dB;
//    } else
//        throw std::runtime_error("Invalid Data");
//}

PullOutput::PullTree_t&PullOutput::getProtonTree(const Fitter::FitParticle& particle)
{
    const auto& det = particle.Particle->Candidate->Detector;

    if(particle.Particle->Type() == ParticleTypeDatabase::Photon) {

        if(det & Detector_t::Type_t::CB) {
            return photons_cb;
        } else if(det & Detector_t::Type_t::TAPS) {
            return photons_taps;
        } else
            throw std::runtime_error("Unexpected detector type in fitter");

    } if(particle.Particle->Type() == ParticleTypeDatabase::Proton) {

        if(det & Detector_t::Type_t::CB) {
            return proton_cb;
        } else if(det & Detector_t::Type_t::TAPS) {
            return proton_taps;
        } else
            throw std::runtime_error("Unexpected detector type in fitter");


    } else
        throw std::runtime_error("Unexpected Particle type in fitter!");
}

PullOutput::PullOutput(HistogramFactory& histfac)
{
    photons_cb.CreateBranches(  histfac.makeTTree("pulls_photon_cb"));
    photons_taps.CreateBranches(histfac.makeTTree("pulls_photon_taps"));
    proton_cb.CreateBranches(   histfac.makeTTree("pulls_proton_cb"));
    proton_taps.CreateBranches( histfac.makeTTree("pulls_proton_taps"));

}

PullOutput::~PullOutput()
{}

void PullOutput::Fill(const std::vector<Fitter::FitParticle>& fitParticles)
{

    for(const auto& p : fitParticles) {

        auto& tree = getProtonTree(p);

        ///@todo replace by initial values
        tree.E     = p.Ek.Value;
        tree.Theta = p.Theta.Value;
        tree.Phi   = p.Phi.Value;

        tree.PullE     = p.Ek.Pull;
        tree.PullTheta = p.Theta.Pull;
        tree.PullPhi   = p.Phi.Pull;

        tree.Tree->Fill();

    }
}
