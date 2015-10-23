#include "ParticleTypeTree.h"

#include <stdexcept>

using namespace std;
using namespace ant;


ParticleTypeTreeDatabase::database_t ParticleTypeTreeDatabase::database = ParticleTypeTreeDatabase::CreateDatabase();
bool ParticleTypeTreeDatabase::is_sorted = false;

void add_Type_2g(ParticleTypeTree& t, const ParticleTypeDatabase::Type& type) {
    auto d = t->CreateDaughter(type);
    d->CreateDaughter(ParticleTypeDatabase::Photon);
    d->CreateDaughter(ParticleTypeDatabase::Photon);
}

void add_Pi0_2g(ParticleTypeTree& t) {
    add_Type_2g(t, ParticleTypeDatabase::Pi0);
}

void add_Pi0_gEpEm(ParticleTypeTree& t) {
    auto d = t->CreateDaughter(ParticleTypeDatabase::Pi0);
    d->CreateDaughter(ParticleTypeDatabase::Photon);
    d->CreateDaughter(ParticleTypeDatabase::ePlus);
    d->CreateDaughter(ParticleTypeDatabase::eMinus);
}

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

ParticleTypeTreeDatabase::database_t ParticleTypeTreeDatabase::CreateDatabase()
{
    database_t database;

    database[Channel::Direct1Pi0_2g] = GetBaseTree();
    add_Pi0_2g(database[Channel::Direct1Pi0_2g]);

    database[Channel::Direct2Pi0_4g] = GetBaseTree();
    add_Pi0_2g(database[Channel::Direct2Pi0_4g]);
    add_Pi0_2g(database[Channel::Direct2Pi0_4g]);

    database[Channel::Direct3Pi0_6g] = GetBaseTree();
    add_Pi0_2g(database[Channel::Direct3Pi0_6g]);
    add_Pi0_2g(database[Channel::Direct3Pi0_6g]);
    add_Pi0_2g(database[Channel::Direct3Pi0_6g]);

    database[Channel::Direct2Pi0_2ggEpEm] = GetBaseTree();
    add_Pi0_2g(database[Channel::Direct2Pi0_2ggEpEm]);
    add_Pi0_gEpEm(database[Channel::Direct2Pi0_2ggEpEm]);

    database[Channel::Direct3Pi0_4ggEpEm] = GetBaseTree();
    add_Pi0_2g(database[Channel::Direct3Pi0_4ggEpEm]);
    add_Pi0_2g(database[Channel::Direct3Pi0_4ggEpEm]);
    add_Pi0_gEpEm(database[Channel::Direct3Pi0_4ggEpEm]);

    database[Channel::DirectPi0Eta_4g] = GetBaseTree();
    add_Pi0_2g(database[Channel::DirectPi0Eta_4g]);
    add_Type_2g(database[Channel::DirectPi0Eta_4g], ParticleTypeDatabase::Eta);

    auto make_Omega_gPseudoscalar_3g = [] (const ParticleTypeDatabase::Type& etapi_type) {
        auto t = GetBaseTree();
        auto omega = t->CreateDaughter(ParticleTypeDatabase::Omega);
        omega->CreateDaughter(ParticleTypeDatabase::Photon);
        add_Type_2g(omega, etapi_type);
        return t;
    };
    database[Channel::Omega_gEta_3g] = make_Omega_gPseudoscalar_3g(ParticleTypeDatabase::Eta);
    database[Channel::Omega_gPi0_3g] = make_Omega_gPseudoscalar_3g(ParticleTypeDatabase::Pi0);

    database[Channel::EtaPrime_2g] = GetBaseTree();
    add_Type_2g(database[Channel::EtaPrime_2g], ParticleTypeDatabase::EtaPrime);


    auto make_EtaPrime_2Pi0Pseudoscalar_6g = [] (const ParticleTypeDatabase::Type& etapi_type) {
        auto t = GetBaseTree();
        auto etap = t->CreateDaughter(ParticleTypeDatabase::EtaPrime);
        add_Type_2g(etap, etapi_type);
        add_Pi0_2g(etap);
        add_Pi0_2g(etap);
        return t;
    };
    database[Channel::EtaPrime_3Pi0_6g] = make_EtaPrime_2Pi0Pseudoscalar_6g(ParticleTypeDatabase::Pi0);
    database[Channel::EtaPrime_2Pi0Eta_6g] = make_EtaPrime_2Pi0Pseudoscalar_6g(ParticleTypeDatabase::Eta);


    auto make_EtaPrime_gOmega_ggPi0_4g = [] () {
        auto t = GetBaseTree();
        auto etap = t->CreateDaughter(ParticleTypeDatabase::EtaPrime);
        etap->CreateDaughter(ParticleTypeDatabase::Photon);
        auto omega = etap->CreateDaughter(ParticleTypeDatabase::Omega);
        omega->CreateDaughter(ParticleTypeDatabase::Photon);
        add_Pi0_2g(omega);
        return t;
    };
    database[Channel::EtaPrime_gOmega_ggPi0_4g] = make_EtaPrime_gOmega_ggPi0_4g();

    // do not sort them here, since references to ParticleTypeDatabase::Type's might not be initialized yet!

    return database;
}

ParticleTypeTree ParticleTypeTreeDatabase::GetBaseTree()
{
    ParticleTypeTree t = Tree<typename ParticleTypeTree::element_type::type>::MakeNode(ParticleTypeDatabase::BeamProton);
    t->CreateDaughter(ParticleTypeDatabase::Proton);
    return t;
}

