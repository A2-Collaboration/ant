#pragma once

#include "tree/TEventData.h"
#include "tree/TCandidate.h"

#include "TH1D.h"


namespace ant {
namespace analysis {
namespace physics {


struct tools
{
    struct protonSelection_t
    {
        TParticlePtr   Proton;
        TParticleList  Photons;
        LorentzVec     PhotonSum;
        LorentzVec     Proton_MM;
        LorentzVec     PhotonBeam;
        double         Copl_pg;
        double         Angle_pMM;
        double         Tagg_E;
        protonSelection_t(const TParticlePtr& proton, const TParticleList& photons,
                          const LorentzVec& photonSum,
                          const LorentzVec& protonMM,
                          const LorentzVec& phtonBeam,
                          double copl_pg,
                          double angle_pMM,
                          double tagg_E):
            Proton(proton),
            Photons(photons),
            PhotonSum(photonSum),
            Proton_MM(protonMM),
            PhotonBeam(phtonBeam),
            Copl_pg(copl_pg),
            Angle_pMM(angle_pMM),
            Tagg_E(tagg_E){}
// neeed?
//        protonSelection_t() = default;

    };


    /// TODO: find out type of TCandidateList::get_iter get rid of template to move to cc....
    template<typename candidateIt>
    static protonSelection_t getProtonSelection(const candidateIt& selectedProton,
                                                const TCandidateList& candidates,
                                                const LorentzVec& photonBeam, double taggE)
    {
        const auto proton = std::make_shared<TParticle>(ParticleTypeDatabase::Proton, selectedProton);
        TParticleList gammas;
        LorentzVec    photonSum;
        for ( auto i_photon : candidates.get_iter())
            if (!(i_photon == selectedProton))
            {
                gammas.emplace_back(std::make_shared<TParticle>(ParticleTypeDatabase::Photon, i_photon));
                photonSum += *gammas.back();
            }
        const auto protonMM = photonBeam + LorentzVec({0, 0, 0}, ParticleTypeDatabase::Proton.Mass())- photonSum;

        return protonSelection_t(proton,
                                 gammas,
                                 photonSum,
                                 protonMM,
                                 photonBeam,
                                 std_ext::radian_to_degree(vec2::Phi_mpi_pi(proton->Phi()-photonSum.Phi() - M_PI )),
                                 std_ext::radian_to_degree(protonMM.Angle(proton->p)),
                                 taggE
                                 );
    }

    template<class T>
    static bool cutOn(const std::string& fillName, const interval<T>& range, const T& val, TH1D* steps)
    {
        auto pass = range.Contains(val);
        if (pass)
        {
            std::string s = std_ext::formatter() << fillName << ": " << range;
            steps->Fill(s.c_str(),1);
        }
        return !pass;
    }

    static double getChargedClusterE(const TClusterList& clusters);
    static double getChargedCandidateE(const TCandidateList& cands);


    //strict stolen from oli:
    static double getCorrVetoEnergy(const TCandidate& photon, const TCandidate& proton);

    static double getPhotonVetoEnergy(const protonSelection_t& sel, const bool strict = false);

};


}
}
}
