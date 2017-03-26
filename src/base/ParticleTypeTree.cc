#include "ParticleTypeTree.h"

#include <stdexcept>

using namespace std;
using namespace ant;

ParticleTypeTreeDatabase::database_t ParticleTypeTreeDatabase::database = ParticleTypeTreeDatabase::CreateDatabase();
bool ParticleTypeTreeDatabase::is_sorted = false;

ParticleTypeTree ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel channel)
{
    // hopefully, at this point the references are initialized
    if(!is_sorted) {
        // sort them all by default operator<
        std::for_each(database.begin(), database.end(), [] (database_t::value_type& item) {
            item.second->Sort();
        });
        is_sorted = true;
    }

    auto it = database.find(channel);
    if(it == database.end())
        throw runtime_error("Did not find ParticleTypeTree for requested channel");
    return it->second;
}

// this little tuple maker is hopefully not a trouble maker :)
template<typename... Types>
std::tuple<typename std::decay<Types>::type...> to(Types&&... types) {
    return std::make_tuple(std::forward<Types>(types)...);
}

ParticleTypeTreeDatabase::database_t ParticleTypeTreeDatabase::CreateDatabase()
{
    database_t database;

    using PTree_t = ParticleTypeTree::element_type; // get underlying tree type from shared_ptr
    auto& gp  = ParticleTypeDatabase::BeamProton;
//    auto& gn  = ParticleTypeDatabase::BeamNeutron; // production off neutron not implemented at all


    auto g    = PTree_t::MakeNode(ParticleTypeDatabase::Photon);
    auto p    = PTree_t::MakeNode(ParticleTypeDatabase::Proton);
    auto n    = PTree_t::MakeNode(ParticleTypeDatabase::Neutron);
    auto eM   = PTree_t::MakeNode(ParticleTypeDatabase::eMinus);
    auto eP   = PTree_t::MakeNode(ParticleTypeDatabase::ePlus);
    auto piP  = PTree_t::MakeNode(ParticleTypeDatabase::PiPlus);
    auto piM  = PTree_t::MakeNode(ParticleTypeDatabase::PiMinus);

    auto pi0   = PTree_t::MakeNode(ParticleTypeDatabase::Pi0);
    auto eta   = PTree_t::MakeNode(ParticleTypeDatabase::Eta);
    auto omega = PTree_t::MakeNode(ParticleTypeDatabase::Omega);
    auto rho   = PTree_t::MakeNode(ParticleTypeDatabase::Rho);
    auto etap  = PTree_t::MakeNode(ParticleTypeDatabase::EtaPrime);
    auto k0s   = PTree_t::MakeNode(ParticleTypeDatabase::K0s);

    auto sigmaPlus = PTree_t::MakeNode(ParticleTypeDatabase::SigmaPlus);
    auto deltaPlus = PTree_t::MakeNode(ParticleTypeDatabase::DeltaPlus);

    // direct multi pi0/eta
    database[Channel::TwoPi0_4g]
            = PTree_t::Make(gp, to(p,
                                   pi0, to(g, g),
                                   pi0, to(g, g)));
    database[Channel::ThreePi0_6g]
            = PTree_t::Make(gp, to(p,
                                   pi0, to(g, g),
                                   pi0, to(g, g),
                                   pi0, to(g, g)));
    database[Channel::TwoPi0_2ggEpEm]
            = PTree_t::Make(gp, to(p,
                                   pi0, to(g, g),
                                   pi0, to(eP, eM, g)));
    database[Channel::ThreePi0_4ggEpEm]
            = PTree_t::Make(gp, to(p,
                                   pi0, to(g, g),
                                   pi0, to(g, g),
                                   pi0, to(eP, eM, g)));;
    database[Channel::Pi0Eta_4g]
            = PTree_t::Make(gp, to(p,
                                   pi0, to(g, g),
                                   eta, to(g, g)));
    database[Channel::Pi0Eta_Pi03Pi0_8g]
            = PTree_t::Make(gp, to(p,
                                   pi0, to(g, g),
                                   eta, to(
                                       pi0, to(g, g),
                                       pi0, to(g, g),
                                       pi0, to(g, g))));
    database[Channel::Pi0Eta_gEpEm2g]
            = PTree_t::Make(gp, to(p,
                                   pi0, to(eP, eM, g),
                                   eta, to(g, g)));
    database[Channel::Pi0Eta_2gPiPi2g]
            = PTree_t::Make(gp, to(p,
                                   pi0, to(g, g),
                                   eta, to(g, g, piP, piM)));
    database[Channel::Pi0PiPi_2gPiPi]
            = PTree_t::Make(gp, to(p,
                                   pi0, to(g, g),
                                   piP,
                                   piM));
    database[Channel::TwoPi0PiPi_4gPiPi]
            = PTree_t::Make(gp, to(p,
                                   pi0, to(g, g),
                                   pi0, to(g, g),
                                   piP,
                                   piM));

    // pi0 decays
    database[Channel::Pi0_2g]
            = PTree_t::Make(gp, to(p,
                                   pi0, to(g, g)));
    database[Channel::Pi0_eeg]
            = PTree_t::Make(gp, to(p,
                                   pi0, to(eP, eM, g)));

    // eta decays
    database[Channel::Eta_2g]
            = PTree_t::Make(gp, to(p,
                                   eta, to(g, g)));
    database[Channel::Eta_eeg]
            = PTree_t::Make(gp, to(p,
                                   eta, to(eP, eM, g)));
    database[Channel::Eta_3Pi0_6g]
            = PTree_t::Make(gp, to(p,
                                   eta, to(
                                       pi0, to(g, g),
                                       pi0, to(g, g),
                                       pi0, to(g, g))));
    database[Channel::Eta_4Pi0_8g]
            = PTree_t::Make(gp, to(p,
                                   eta, to(
                                       pi0, to(g, g),
                                       pi0, to(g, g),
                                       pi0, to(g, g),
                                       pi0, to(g, g))));

    // omega decays
    database[Channel::Omega_Pi0PiPPiM_2g]
            = PTree_t::Make(gp, to(p,
                                   omega, to(
                                       pi0, to(g, g),
                                       piP,
                                       piM)));
    database[Channel::Omega_gEta_3g]
            = PTree_t::Make(gp, to(p,
                                   omega, to(
                                       g,
                                       eta, to(g, g))));
    database[Channel::Omega_gPi0_3g]
            = PTree_t::Make(gp, to(p,
                                   omega, to(
                                       g,
                                       pi0, to(g, g))));

    // etaprime decays
    database[Channel::EtaPrime_2g]
            = PTree_t::Make(gp, to(p,
                                   etap, to(g, g)));
    database[Channel::EtaPrime_3Pi0_6g]
            = PTree_t::Make(gp, to(p,
                                   etap, to(
                                       pi0, to(g, g),
                                       pi0, to(g, g),
                                       pi0, to(g, g))));
    database[Channel::EtaPrime_2Pi0Eta_6g]
            = PTree_t::Make(gp, to(p,
                                   etap, to(
                                       pi0, to(g, g),
                                       pi0, to(g, g),
                                       eta, to(g, g))));
    database[Channel::EtaPrime_EtaPiPPiM_2gPiPPiM]
            = PTree_t::Make(gp, to(p,
                                   etap, to(
                                       eta, to(g, g),
                                       piP,
                                       piM)));
    database[Channel::EtaPrime_gOmega_ggPi0_4g]
            = PTree_t::Make(gp, to(p,
                                   etap, to(
                                       g,
                                       omega, to(
                                           g,
                                           pi0, to(g, g)))));
    database[Channel::EtaPrime_eeg]
            = PTree_t::Make(gp, to(p,
                                   etap, to(eP, eM, g)));

    database[Channel::EtaPrime_gRho_gPiPi]
            = PTree_t::Make(gp, to(p,
                                   etap, to(
                                       g,
                                       rho, to(piP, piM))));

    // rho decays
    database[Channel::Rho_PiPi]
            = PTree_t::Make(gp, to(p,
                                   rho, to(piP, piM)));

    // nucleon resonances
    database[Channel::SigmaPlusK0s_6g]
            = PTree_t::Make(gp, to(
                                k0s, to(
                                    pi0, to(g, g),
                                    pi0, to(g, g)),
                                sigmaPlus, to(
                                    p,
                                    pi0, to(g, g))));
    database[Channel::gp_DeltaPlus2Pi0_3Pi0_6g]
            = PTree_t::Make(gp, to(
                                pi0, to(g, g),
                                pi0, to(g, g),
                                deltaPlus, to(
                                    p,
                                    pi0, to(g, g))));

    // production channels
    database[Channel::gp_pPi0]
            = PTree_t::Make(gp, to(p, pi0));
    database[Channel::gp_pPi0Pi0]
            = PTree_t::Make(gp, to(p, pi0, pi0));
    database[Channel::gp_p3Pi0]
            = PTree_t::Make(gp, to(p, pi0, pi0, pi0));
    database[Channel::gp_pPiPPiMPi0]
            = PTree_t::Make(gp, to(p, piP, piM, pi0));
    database[Channel::gp_pPiPPiMPi0Pi0]
            = PTree_t::Make(gp, to(p, piP, piM, pi0, pi0));
    database[Channel::gp_pEta]
            = PTree_t::Make(gp, to(p, eta));
    database[Channel::gp_pEtaPi0]
            = PTree_t::Make(gp, to(p, eta, pi0));
    database[Channel::gp_pOmega]
            = PTree_t::Make(gp, to(p, omega));
    database[Channel::gp_pEtaPrime]
            = PTree_t::Make(gp, to(p, etap));
    database[Channel::gp_pRho]
            = PTree_t::Make(gp, to(p, rho));
    database[Channel::gp_SigmaPlusK0S]
            = PTree_t::Make(gp, to(sigmaPlus, k0s));
    database[Channel::gp_nPiP]
            = PTree_t::Make(gp, to(n, piP));
    database[Channel::gp_pg]
            = PTree_t::Make(gp, to(p, g));

    // do not sort them here, since references to ParticleTypeDatabase::Type's might not be initialized yet!

    return database;
}


