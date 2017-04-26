#include "tools.h"

#include "base/std_ext/math.h"

#include <numeric>

using namespace std;
using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::physics;

const double strictVetoCutAngle  = 15.0;

template<typename candidateIt>
tools::protonSelection_t getProtonSelection(const candidateIt& selectedProton,
                                            const TCandidateList& candidates,
                                            const LorentzVec& photonBeam, double taggE)
{
    auto protonVeto  = 0.;
    auto photonVeto  = 0.;
    const auto proton = std::make_shared<TParticle>(ParticleTypeDatabase::Proton, selectedProton);
    TParticleList gammas;
    LorentzVec    photonSum;
    TCandidatePtrList photonCands;

    for ( auto i_photon : candidates.get_iter())
    {
        if (!(i_photon == selectedProton))
        {
            photonVeto += i_photon->VetoEnergy;
            gammas.emplace_back(std::make_shared<TParticle>(ParticleTypeDatabase::Photon, i_photon));
            photonSum += *gammas.back();
            photonCands.emplace_back(i_photon);
        }
    }
    protonVeto = selectedProton->VetoEnergy;


    const auto protonMM = photonBeam + LorentzVec({0, 0, 0}, ParticleTypeDatabase::Proton.Mass())- photonSum;

    return tools::protonSelection_t(
                proton, gammas,
                selectedProton, photonCands,
                photonSum,
                protonMM,
                photonBeam,
                std_ext::radian_to_degree(vec2::Phi_mpi_pi(proton->Phi()-photonSum.Phi() - M_PI )),
                std_ext::radian_to_degree(protonMM.Angle(proton->p)),
                taggE,
                protonVeto,
                photonVeto
                );
}

std::vector<tools::protonSelection_t> tools::makeProtonSelections(const TCandidateList& candidates, const LorentzVec& photonBeam, double taggE, const IntervalD& imMMprotonCut)
{
    std::vector<tools::protonSelection_t> psels;
    for (auto i_proton: candidates.get_iter())
    {
        const auto select = getProtonSelection(i_proton,candidates,
                                               photonBeam, taggE);
        if (imMMprotonCut.Contains(select.Proton_MM.M()))
            psels.emplace_back(select);
    }
    return psels;
}

TCandidatePtrList tools::getNeutral(const TEventData& data, const double threshold)
{

    TCandidatePtrList neutrals;
    const auto& cds = data.Candidates;
    for (const auto& c: cds.get_iter())
    {
        if (c->VetoEnergy > threshold)
            continue;
        neutrals.emplace_back(c);
    }
    return neutrals;
}

double tools::getChargedClusterE(const TClusterList& clusters)
{
    auto accE = 0.0;
    for(const auto& cl : clusters)
    {
        if(cl.DetectorType == Detector_t::Type_t::PID ||
           cl.DetectorType == Detector_t::Type_t::TAPSVeto)
            accE += cl.Energy;
    }
    return accE;
}

double tools::getCandidateVetoE(const TCandidateList& cands)
{
    if (cands.size() == 0) return 0;

    return accumulate(cands.begin(),cands.end(),
                      cands[0].VetoEnergy,
            [] (double acc, const TCandidate& ca)
    {
        return acc + ca.VetoEnergy;
    });
}

double tools::getCorrVetoEnergy(const TCandidate& photon, const TCandidate& proton)
{
    if (fabs(vec2::Phi_mpi_pi(proton.Phi - photon.Phi - M_PI)) > degree_to_radian(strictVetoCutAngle))
        return 0;
    return photon.VetoEnergy;
}


double tools::getPhotonVetoEnergy(const tools::protonSelection_t& sel, const bool strict)
{
    double acc = 0;
    for (const auto& ph: sel.Photons)
    {
        acc += strict ? getCorrVetoEnergy(*ph->Candidate,*sel.Proton->Candidate) : ph->Candidate->VetoEnergy;
    }
    return acc;
}

