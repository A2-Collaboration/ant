#include "ParticleTypeTree.h"

#include <stdexcept>

using namespace std;
using namespace ant;

ParticleTypeTree GetBaseTree();
ParticleTypeTree GetProductionTree(const ParticleTypeDatabase::Type& pseudoBeam,
                                   std::vector<std::reference_wrapper<const ParticleTypeDatabase::Type>> products);

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

void add_Type_Dalitz(ParticleTypeTree& t, const ParticleTypeDatabase::Type& type) {
    auto d = t->CreateDaughter(type);
    d->CreateDaughter(ParticleTypeDatabase::Photon);
    d->CreateDaughter(ParticleTypeDatabase::ePlus);
    d->CreateDaughter(ParticleTypeDatabase::eMinus);
}

void add_Pi0_gEpEm(ParticleTypeTree& t) {
    add_Type_Dalitz(t, ParticleTypeDatabase::Pi0);
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

    database[Channel::Pi0_2g] = GetBaseTree();
    add_Pi0_2g(database[Channel::Pi0_2g]);

    database[Channel::Pi0_eeg] = GetBaseTree();
    add_Type_Dalitz(database[Channel::Pi0_eeg], ParticleTypeDatabase::Pi0);

    database[Channel::TwoPi0_4g] = GetBaseTree();
    add_Pi0_2g(database[Channel::TwoPi0_4g]);
    add_Pi0_2g(database[Channel::TwoPi0_4g]);

    database[Channel::ThreePi0_6g] = GetBaseTree();
    add_Pi0_2g(database[Channel::ThreePi0_6g]);
    add_Pi0_2g(database[Channel::ThreePi0_6g]);
    add_Pi0_2g(database[Channel::ThreePi0_6g]);

    database[Channel::TwoPi0_2ggEpEm] = GetBaseTree();
    add_Pi0_2g(database[Channel::TwoPi0_2ggEpEm]);
    add_Pi0_gEpEm(database[Channel::TwoPi0_2ggEpEm]);

    database[Channel::ThreePi0_4ggEpEm] = GetBaseTree();
    add_Pi0_2g(database[Channel::ThreePi0_4ggEpEm]);
    add_Pi0_2g(database[Channel::ThreePi0_4ggEpEm]);
    add_Pi0_gEpEm(database[Channel::ThreePi0_4ggEpEm]);

    database[Channel::Pi0Eta_4g] = GetBaseTree();
    add_Pi0_2g(database[Channel::Pi0Eta_4g]);
    add_Type_2g(database[Channel::Pi0Eta_4g], ParticleTypeDatabase::Eta);

    {
        database[Channel::Pi0Eta_Pi03Pi0_8g] = GetBaseTree();
        auto& base = database[Channel::Pi0Eta_Pi03Pi0_8g];
        add_Pi0_2g(base);
        auto& eta = base->CreateDaughter(ParticleTypeDatabase::Eta);
        add_Pi0_2g(eta);
        add_Pi0_2g(eta);
        add_Pi0_2g(eta);
    }
    {
        database[Channel::Pi0Eta_gEpEm2g] = GetBaseTree();
        auto& base = database[Channel::Pi0Eta_gEpEm2g];
        add_Pi0_gEpEm(base);
        add_Type_2g(base, ParticleTypeDatabase::Eta);
    }
    {
        database[Channel::Pi0Eta_2gPiPi2g] = GetBaseTree();
        auto& base = database[Channel::Pi0Eta_2gPiPi2g];
        add_Pi0_2g(base);
        auto& eta = base->CreateDaughter(ParticleTypeDatabase::Eta);
        eta->CreateDaughter(ParticleTypeDatabase::Photon);
        eta->CreateDaughter(ParticleTypeDatabase::Photon);
        eta->CreateDaughter(ParticleTypeDatabase::PiPlus);
        eta->CreateDaughter(ParticleTypeDatabase::PiMinus);
    }
    {
        database[Channel::Pi0PiPi_2gPiPi] = GetBaseTree();
        auto& base = database[Channel::Pi0PiPi_2gPiPi];
        add_Pi0_2g(base);
        base->CreateDaughter(ParticleTypeDatabase::PiPlus);
        base->CreateDaughter(ParticleTypeDatabase::PiMinus);
    }
    {
        database[Channel::TwoPi0PiPi_4gPiPi] = GetBaseTree();
        auto& base = database[Channel::TwoPi0PiPi_4gPiPi];
        add_Pi0_2g(base);
        add_Pi0_2g(base);
        base->CreateDaughter(ParticleTypeDatabase::PiPlus);
        base->CreateDaughter(ParticleTypeDatabase::PiMinus);
    }

    database[Channel::Eta_2g] = GetBaseTree();
    add_Type_2g(database[Channel::Eta_2g],ParticleTypeDatabase::Eta);

    database[Channel::Eta_eeg] = GetBaseTree();
    add_Type_Dalitz(database[Channel::Eta_eeg], ParticleTypeDatabase::Eta);

    {
        database[Channel::Omega_Pi0PiPPiM_2g] = GetBaseTree();
        auto& base = database[Channel::Omega_Pi0PiPPiM_2g];
        auto omega = base->CreateDaughter(ParticleTypeDatabase::Omega);
        add_Pi0_2g(omega);
        omega->CreateDaughter(ParticleTypeDatabase::PiPlus);
        omega->CreateDaughter(ParticleTypeDatabase::PiMinus);
    }


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


    auto make_Eta_3Pi0_ng = [] (const unsigned n)
    {
        auto etaTree = GetBaseTree();
        auto eta = etaTree->CreateDaughter(ParticleTypeDatabase::Eta);
        for ( auto i = 0u ; i < n ; ++i)
            add_Pi0_2g(eta);
        return etaTree;
    };
    database[Channel::Eta_3Pi0_6g] = make_Eta_3Pi0_ng(3);
    database[Channel::Eta_4Pi0_8g] = make_Eta_3Pi0_ng(4);



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

    auto make_EtaPrime_EtaPiPPiM_2gPiPPiM = [] {
        auto t = GetBaseTree();
        auto etap = t->CreateDaughter(ParticleTypeDatabase::EtaPrime);
        add_Type_2g(etap, ParticleTypeDatabase::Eta);
        etap->CreateDaughter(ParticleTypeDatabase::PiMinus);
        etap->CreateDaughter(ParticleTypeDatabase::PiPlus);
        return t;
    };
    database[Channel::EtaPrime_EtaPiPPiM_2gPiPPiM] = make_EtaPrime_EtaPiPPiM_2gPiPPiM();

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


    database[Channel::EtaPrime_eeg] = GetBaseTree();
    add_Type_Dalitz(database[Channel::EtaPrime_eeg], ParticleTypeDatabase::EtaPrime);

    auto make_EtaPrime_gRho_gPiPi = [] () {
        auto t = GetBaseTree();
        auto etap = t->CreateDaughter(ParticleTypeDatabase::EtaPrime);
        etap->CreateDaughter(ParticleTypeDatabase::Photon);
        auto rho = etap->CreateDaughter(ParticleTypeDatabase::Rho);
        rho->CreateDaughter(ParticleTypeDatabase::PiPlus);
        rho->CreateDaughter(ParticleTypeDatabase::PiMinus);
        return t;
    };
    database[Channel::EtaPrime_gRho_gPiPi] = make_EtaPrime_gRho_gPiPi();


    auto make_Rho_PiPi = [] () {
        auto t = GetBaseTree();
        auto rho = t->CreateDaughter(ParticleTypeDatabase::Rho);
        rho->CreateDaughter(ParticleTypeDatabase::PiPlus);
        rho->CreateDaughter(ParticleTypeDatabase::PiMinus);
        return t;
    };
    database[Channel::Rho_PiPi] = make_Rho_PiPi();

    auto make_SigmaPlusK0s_6g = []
    {
        auto sigmaTree = Tree<typename ParticleTypeTree::element_type::type>::MakeNode(ParticleTypeDatabase::BeamProton);
        auto SigmaPlus = sigmaTree->CreateDaughter(ParticleTypeDatabase::SigmaPlus);
        SigmaPlus->CreateDaughter(ParticleTypeDatabase::Proton);
        add_Pi0_2g(SigmaPlus);
        auto k0Short = sigmaTree->CreateDaughter(ParticleTypeDatabase::K0s);
        add_Pi0_2g(k0Short);
        add_Pi0_2g(k0Short);
        return sigmaTree;
    };
    database[Channel::SigmaPlusK0s_6g] = make_SigmaPlusK0s_6g();

    database[Channel::gp_DeltaPlus2Pi0_3Pi0_6g] = []
    {
        auto t = GetBaseTree();
        auto DeltaPlus = t->CreateDaughter(ParticleTypeDatabase::DeltaPlus);
        DeltaPlus->CreateDaughter(ParticleTypeDatabase::Proton);
        add_Pi0_2g(DeltaPlus);
        add_Pi0_2g(t);
        add_Pi0_2g(t);
        return t;
    }();

    database[Channel::gp_pPi0]          = GetProductionTree(ParticleTypeDatabase::BeamProton,
                                                   {ParticleTypeDatabase::Proton,
                                                    ParticleTypeDatabase::Pi0});
    database[Channel::gp_pPi0Pi0]       = GetProductionTree(ParticleTypeDatabase::BeamProton,
                                                   {ParticleTypeDatabase::Proton,
                                                    ParticleTypeDatabase::Pi0,
                                                    ParticleTypeDatabase::Pi0});
    database[Channel::gp_p3Pi0]       = GetProductionTree(ParticleTypeDatabase::BeamProton,
                                                   {ParticleTypeDatabase::Proton,
                                                    ParticleTypeDatabase::Pi0,
                                                    ParticleTypeDatabase::Pi0,
                                                    ParticleTypeDatabase::Pi0});
    database[Channel::gp_pPiPPiMPi0]    = GetProductionTree(ParticleTypeDatabase::BeamProton,
                                                   {ParticleTypeDatabase::Proton,
                                                    ParticleTypeDatabase::PiPlus,
                                                    ParticleTypeDatabase::PiMinus,
                                                    ParticleTypeDatabase::Pi0});
    database[Channel::gp_pPiPPiMPi0Pi0] = GetProductionTree(ParticleTypeDatabase::BeamProton,
                                                   {ParticleTypeDatabase::Proton,
                                                    ParticleTypeDatabase::PiPlus,
                                                    ParticleTypeDatabase::PiMinus,
                                                    ParticleTypeDatabase::Pi0,
                                                    ParticleTypeDatabase::Pi0});
    database[Channel::gp_pg]          = GetProductionTree(ParticleTypeDatabase::BeamProton,
                                                   {ParticleTypeDatabase::Proton,
                                                    ParticleTypeDatabase::Photon});
    database[Channel::gp_nPiP]          = GetProductionTree(ParticleTypeDatabase::BeamProton,
                                                   {ParticleTypeDatabase::Neutron,
                                                    ParticleTypeDatabase::PiPlus});
    database[Channel::gp_pEta]          = GetProductionTree(ParticleTypeDatabase::BeamProton,
                                                   {ParticleTypeDatabase::Proton,
                                                    ParticleTypeDatabase::Eta});
    database[Channel::gp_pEtaPi0]       = GetProductionTree(ParticleTypeDatabase::BeamProton,
                                                   {ParticleTypeDatabase::Proton,
                                                    ParticleTypeDatabase::Eta,
                                                    ParticleTypeDatabase::Pi0});
    database[Channel::gp_pEtaPrime]     = GetProductionTree(ParticleTypeDatabase::BeamProton,
                                                   {ParticleTypeDatabase::Proton,
                                                    ParticleTypeDatabase::EtaPrime});
    database[Channel::gp_pOmega]        = GetProductionTree(ParticleTypeDatabase::BeamProton,
                                                   {ParticleTypeDatabase::Proton,
                                                    ParticleTypeDatabase::Omega});
    database[Channel::gp_pRho]          = GetProductionTree(ParticleTypeDatabase::BeamProton,
                                                   {ParticleTypeDatabase::Proton,
                                                    ParticleTypeDatabase::Rho});
    database[Channel::gp_SigmaPlusK0S]  = GetProductionTree(ParticleTypeDatabase::BeamProton,
                                                   {ParticleTypeDatabase::SigmaPlus,
                                                    ParticleTypeDatabase::K0s});


    // do not sort them here, since references to ParticleTypeDatabase::Type's might not be initialized yet!

    return database;
}

ParticleTypeTree GetBaseTree()
{
    ParticleTypeTree t = Tree<typename ParticleTypeTree::element_type::type>::MakeNode(ParticleTypeDatabase::BeamProton);
    t->CreateDaughter(ParticleTypeDatabase::Proton);
    return t;
}

ParticleTypeTree GetProductionTree(const ParticleTypeDatabase::Type& pseudoBeam,
                                   std::vector<std::reference_wrapper<const ParticleTypeDatabase::Type>> products)
{
    auto t = Tree<typename ParticleTypeTree::element_type::type>::MakeNode(pseudoBeam);
    for (auto& p: products)
    {
        t->CreateDaughter(p);
    }
    return t;
}

