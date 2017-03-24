#pragma once

#include "ParticleType.h"
#include "Tree.h"

#include <map>

namespace ant {

using ParticleTypeTree = std::shared_ptr<Tree<const ParticleTypeDatabase::Type&>>;

class ParticleTypeTreeDatabase {
public:
    enum class Channel {
        ThreePi0_6g,
        TwoPi0_4g,
        Pi0_2g,
        Pi0_eeg,
        TwoPi0_2ggEpEm,
        ThreePi0_4ggEpEm,
        Pi0Eta_4g,
        Pi0Eta_2gPiPi2g,
        Pi0Eta_Pi03Pi0_8g,
        Pi0Eta_gEpEm2g,
        Pi0PiPi_2gPiPi,
        TwoPi0PiPi_4gPiPi,
        Eta_2g,
        Eta_eeg,
        Eta_3Pi0_6g,
        Eta_4Pi0_8g,
        Omega_gEta_3g,
        Omega_gPi0_3g,
        Omega_Pi0PiPPiM_2g,
        EtaPrime_2g,
        EtaPrime_3Pi0_6g,
        EtaPrime_2Pi0Eta_6g,
        EtaPrime_gOmega_ggPi0_4g,
        EtaPrime_eeg,
        EtaPrime_gRho_gPiPi,
        EtaPrime_EtaPiPPiM_2gPiPPiM,
        Rho_PiPi,
        SigmaPlusK0s_6g,
        gp_DeltaPlus2Pi0_3Pi0_6g,
        // production trees used by MC cocktail
        gp_pPi0,
        gp_pPi0Pi0,
        gp_p3Pi0,
        gp_pPiPPiMPi0,
        gp_pPiPPiMPi0Pi0,
        gp_pg,
        gp_nPiP,
        gp_pEta,
        gp_pEtaPi0,
        gp_pEtaPrime,
        gp_pOmega,
        gp_pRho,
        gp_SigmaPlusK0S
    };

    static ParticleTypeTree Get(Channel channel);

protected:
    using database_t = std::map<Channel, ParticleTypeTree>;
    static database_t database;
    static bool is_sorted;

    static database_t CreateDatabase();

public:

    class const_iterator : public database_t::const_iterator {
    public:
        const_iterator(const database_t::const_iterator& i) : database_t::const_iterator(i)  {}
        Channel operator*() const { return database_t::const_iterator::operator*().first; }
    };

    /*
     * Use begin()/end() in for-ranged loop to iterate all Channel enums as follows:
     *
     * for(auto channel : ParticleTypeTreeDatabase()) {
     *   // type is reference to ParticleTypeTreeDatabase::Channel
     *   auto ptree = ParticleTypeTreeDatabase::Get(channel);
     * }
     */
    static const_iterator begin() { return const_iterator(database.begin()); }
    static const_iterator end()   { return const_iterator(database.end()); }


};

}
