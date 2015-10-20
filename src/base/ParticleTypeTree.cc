#include "ParticleTypeTree.h"


using namespace std;
using namespace ant;


ParticleTypeTreeDatabase::database_t ParticleTypeTreeDatabase::database = ParticleTypeTreeDatabase::CreateDatabase();

void add_Type_2g(ParticleTypeTree& t, const ParticleTypeDatabase::Type& type) {
    auto d = t->CreateDaughter(type);
    d->CreateDaughter(ParticleTypeDatabase::Photon);
    d->CreateDaughter(ParticleTypeDatabase::Photon);
}

void add_Pi0_2g(ParticleTypeTree& t) {
    add_Type_2g(t, ParticleTypeDatabase::Pi0);
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


    return database;
}

ParticleTypeTree ParticleTypeTreeDatabase::GetBaseTree()
{
    ParticleTypeTree t = Tree<typename ParticleTypeTree::element_type::type>::MakeNode(ParticleTypeDatabase::BeamProton);
    t->CreateDaughter(ParticleTypeDatabase::Proton);
    return t;
}

