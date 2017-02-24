#include "base/ParticleType.h"
#include <iostream>

#include "base/std_ext/math.h"

using namespace std;
using namespace ant;

namespace ant {
ostream& operator<<(ostream &stream, const ParticleTypeDatabase::Type& particle_type)
{
    stream << "ParticleType " << particle_type.Name() << "(";
    stream << "Mass=" << particle_type.Mass();
    stream << (particle_type.Charged() ? "MeV, Charged)" : "MeV, Neutral)");
    return stream;
}
}

unsigned ParticleTypeDatabase::Type::NextUID;
ParticleTypeDatabase::Particles_t ParticleTypeDatabase::types;

const ParticleTypeDatabase::Type ParticleTypeDatabase::Nucleon("Nucleon",          "Nucleon",     "unknown",   std_ext::NaN, false);
const ParticleTypeDatabase::Type ParticleTypeDatabase::Proton("Proton",            "p",           "p",         938.272046, true, &ParticleTypeDatabase::Nucleon);
const ParticleTypeDatabase::Type ParticleTypeDatabase::Neutron("Neutron",          "n",           "n",         939.565378, false, &ParticleTypeDatabase::Nucleon);
const ParticleTypeDatabase::Type ParticleTypeDatabase::SigmaPlus("SigmaPlus",      "#Sigma^{+}",  "Sigma+",    1189.37,    true);
const ParticleTypeDatabase::Type ParticleTypeDatabase::Photon("Photon",            "#gamma",      "g",         0.0,        false);

const ParticleTypeDatabase::Type ParticleTypeDatabase::Pi0("Pi0",                  "#pi^{0}",     "pi0",       134.9766,  false);
const ParticleTypeDatabase::Type ParticleTypeDatabase::PiCharged("PiCharged",      "#pi^{#pm}",   "unknown",   139.57018, true);
const ParticleTypeDatabase::Type ParticleTypeDatabase::PiPlus("PiPlus",            "#pi^{+}",     "pi+"    ,   139.57018, true, &ParticleTypeDatabase::PiCharged);
const ParticleTypeDatabase::Type ParticleTypeDatabase::PiMinus("PiMinus",          "#pi^{-}",     "pi-",       139.57018, true, &ParticleTypeDatabase::PiCharged);

const ParticleTypeDatabase::Type ParticleTypeDatabase::K0s("K0s",                  "#k^{0}_s",    "K0S",       497.614,   false);


const ParticleTypeDatabase::Type ParticleTypeDatabase::eCharged("eCharged",        "e^{#pm}",     "unknown",   0.510998928, true);
const ParticleTypeDatabase::Type ParticleTypeDatabase::ePlus("Positron",           "e^{+}",       "e+",        0.510998928, true, &ParticleTypeDatabase::eCharged);
const ParticleTypeDatabase::Type ParticleTypeDatabase::eMinus("Electron",          "e^{-}",       "e-",        0.510998928, true, &ParticleTypeDatabase::eCharged);

const ParticleTypeDatabase::Type ParticleTypeDatabase::MuCharged("MuCharged",      "#mu^{#pm}",   "unknown",   105.658389, true);
const ParticleTypeDatabase::Type ParticleTypeDatabase::MuPlus("MuPlus",            "#mu^{+}",     "mu+",       105.658389, true, &ParticleTypeDatabase::MuCharged);
const ParticleTypeDatabase::Type ParticleTypeDatabase::MuMinus("MuMinus",          "#mu^{-}",     "mu-",       105.658389, true, &ParticleTypeDatabase::MuCharged);

const ParticleTypeDatabase::Type ParticleTypeDatabase::Eta("Eta",                  "#eta",        "eta",       547.853, false);
const ParticleTypeDatabase::Type ParticleTypeDatabase::Omega("Omega",              "#omega",      "omega",     782.65, false);
const ParticleTypeDatabase::Type ParticleTypeDatabase::EtaPrime("EtaPrime",        "#eta'",       "eta'",      957.78, false);
const ParticleTypeDatabase::Type ParticleTypeDatabase::Rho("Rho",                  "#rho'",       "rho0",      775.26, false);

const ParticleTypeDatabase::Type ParticleTypeDatabase::BeamTarget("BeamTarget",    "BeamTarget",  "unknwon",   std_ext::NaN, false);
const ParticleTypeDatabase::Type ParticleTypeDatabase::BeamProton("BeamProton",    "(#gamma p)",  "g p",  938.272046, true,  &ParticleTypeDatabase::BeamTarget);
const ParticleTypeDatabase::Type ParticleTypeDatabase::BeamNeutron("BeamNeutron",  "(#gamma n)",  "g n",  939.565378, false, &ParticleTypeDatabase::BeamTarget);

ParticleTypeDatabase::Type::Type(const string &_name,
                                 const string &_print_name, const string& _pluto_name,
                                 const mev_t &_mass, const bool &_charged, const ParticleTypeDatabase::Type *_sametype):
    UID(NextUID++),
    name(_name),
    print_name(_print_name),
    pluto_name(_pluto_name),
    mass(_mass),
    charged(_charged),
    sametype(_sametype)
{
    types.emplace(UID, *this);
}

void ParticleTypeDatabase::Print()
{
    for(auto& p : types) {
        cout << p.second << endl;
    }
}

const ParticleTypeDatabase::Type *ParticleTypeDatabase::GetTypeOfPlutoID(index_t pid)
{
    PlutoIDMap_t::const_iterator entry = pluto_id_map.find(pid);
    if(entry == pluto_id_map.end()) {
        return nullptr;
    }
    return entry->second;
}

mev_t ParticleTypeDatabase::CalculatePhotoproductionThreshold(mev_t m_sum, const Type& target)
{
    // assumes target with non-zero rest mass
    return m_sum*(m_sum + 2*target.Mass())/(2*target.Mass());
}

const ParticleTypeDatabase::TypeList_t ParticleTypeDatabase::detectables = { &ParticleTypeDatabase::Photon,
                                                                             &ParticleTypeDatabase::Proton,
                                                                             &ParticleTypeDatabase::PiCharged,
                                                                             &ParticleTypeDatabase::eCharged
                                                                           };

const ParticleTypeDatabase::TypeList_t ParticleTypeDatabase::mc_finalstate = { &ParticleTypeDatabase::Photon,
                                                                               &ParticleTypeDatabase::Proton,
                                                                               &ParticleTypeDatabase::PiMinus,
                                                                               &ParticleTypeDatabase::PiPlus,
                                                                               &ParticleTypeDatabase::eMinus,
                                                                               &ParticleTypeDatabase::ePlus,
                                                                               &ParticleTypeDatabase::MuMinus,
                                                                               &ParticleTypeDatabase::MuPlus,
                                                                             };

const ParticleTypeDatabase::TypeList_t ParticleTypeDatabase::neutral_mesons = { &ParticleTypeDatabase::Pi0,
                                                                                &ParticleTypeDatabase::Eta,
                                                                                &ParticleTypeDatabase::Omega,
                                                                                &ParticleTypeDatabase::EtaPrime
                                                                              };


ParticleTypeDatabase::PlutoIDMap_t ParticleTypeDatabase::pluto_id_map = [] () {
    ParticleTypeDatabase::PlutoIDMap_t m;
    m[1]  = &ParticleTypeDatabase::Photon;
    m[2]  = &ParticleTypeDatabase::ePlus;
    m[3]  = &ParticleTypeDatabase::eMinus;
    m[5]  = &ParticleTypeDatabase::MuPlus;
    m[6]  = &ParticleTypeDatabase::MuMinus;
    m[7]  = &ParticleTypeDatabase::Pi0;
    m[8]  = &ParticleTypeDatabase::PiPlus;
    m[9]  = &ParticleTypeDatabase::PiMinus;
    m[13] = &ParticleTypeDatabase::Neutron;
    m[14] = &ParticleTypeDatabase::Proton;
    m[16] = &ParticleTypeDatabase::K0s;
    m[17] = &ParticleTypeDatabase::Eta;
    m[19] = &ParticleTypeDatabase::SigmaPlus;
    m[41] = &ParticleTypeDatabase::Rho;
    m[53] = &ParticleTypeDatabase::EtaPrime;
    m[52] = &ParticleTypeDatabase::Omega;
    m[14001] = &ParticleTypeDatabase::BeamProton;
    return m;
}();
