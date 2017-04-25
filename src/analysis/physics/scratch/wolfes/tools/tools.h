#pragma once

#include "tree/TEventData.h"
#include "tree/TCandidate.h"

#include "TH1D.h"
#include "TLorentzVector.h"


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

    static std::vector<protonSelection_t> makeProtonSelections(
            const TCandidateList& candidates,
            const LorentzVec& photonBeam,
            double taggE,
            const IntervalD& imMMprotonCut = {-std_ext::inf,std_ext::inf}
            );

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

    static std::vector<TLorentzVector> MakeTLorenz(const TParticleList& particles)
    {
        std::vector<TLorentzVector> lg(particles.size());
        std::transform(particles.begin(),particles.end(),lg.begin(),
                       [](const TParticlePtr& ph){return TLorentzVector(*ph);});
        return lg;
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
