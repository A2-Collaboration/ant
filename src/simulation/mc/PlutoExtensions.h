#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#pragma GCC diagnostic ignored "-Wvla"
#pragma GCC diagnostic ignored "-Wignored-qualifiers"
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic ignored "-Wwrite-strings"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include "PStaticData.h"
#pragma GCC diagnostic pop

#include "iostream"

namespace ant
{
namespace simulation
{
namespace mc
{


/**
 * @brief UpdatePluteDataBase:  Place your any additions to Pluto here.
 *                              This function should be called at the
 *        very beginnig when including any Pluto-Objects using the
 *        Pluto-Database.
 */
void UpdatePluteDataBase()
{
    PStaticData* sdata = makeStaticData();


    std::cout << "=============================================" << endl;
    std::cout << "omega before change:" << endl;
    listParticle("w");

    // additional omega decays (from PDG Booklet 2014)
    sdata->AddDecay("w --> eta + g",        "w", "eta,g",           4.6E-4);
    sdata->AddDecay("w --> g + g + g",      "w", "g,g,g",           1.9E-4);    // upper limit
    sdata->AddDecay("w --> pi0 e+ e-",      "w", "pi0,e+,e-",       7.7E-4);
    sdata->AddDecay("w --> pi0 mu+ mu-",    "w", "pi0,mu+,mu-",     1.3E-4);
    sdata->AddDecay("w --> pi+ pi- pi0 pi0","w", "pi+,pi-,pi0,pi0", 2.0E-4);    // upper limit
    sdata->AddDecay("w --> pi+ pi- pi+ pi-","w", "pi+,pi-,pi+,pi-", 1.0E-3);    // upper limit
    sdata->AddDecay("w --> pi0 pi0 g",      "w", "pi0,pi0,g",       6.6E-5);
    sdata->AddDecay("w --> eta pi0 g",      "w", "eta,pi0,g",       3.3E-5);    // upper limit

    // Charge conjucation violating modes
    sdata->AddDecay("w --> eta pi0",        "w", "eta,pi0",         2.1E-4);    // upper limit
    sdata->AddDecay("w --> pi0 pi0",        "w", "pi0,pi0",         2.1E-4);    // upper limit
    sdata->AddDecay("w --> pi0 pi0 pi0",    "w", "pi0,pi0,pi0",     2.3E-4);    // upper limit


    std::cout << "=============================================" << endl;
    std::cout << "omega after change:" << endl;
    listParticle("w");

    std::cout << "=============================================" << endl;
    std::cout << "eta' before change:" << endl;
    listParticle("eta'");
    // newer Decay BRs for eta' (PDG particle listings)
    sdata->SetDecayBR("eta'", "pi+,pi-,eta",    0.429,    1);   // Mode flag:   ( see Pluto-source)
    sdata->SetDecayBR("eta'", "rho0,g",         0.291,    1);   // 0: Add the new b.r. to the existing ones + re-weighting
    sdata->SetDecayBR("eta'", "eta,pi0,pi0",    0.222,    1);   // 1: No re-weighting (in this case br must be <1.)
    sdata->SetDecayBR("eta'", "w,g",            0.0275,   1);
    sdata->SetDecayBR("eta'", "g,g",            0.0220,   1);
    sdata->SetDecayBR("eta'", "pi0,pi0,pi0",    0.00214,  1);
    sdata->SetDecayBR("eta'", "dimuon,g",       0.000108, 1);

    // new Decays from PDG
    sdata->AddDecay("eta' --> pi+ pi- mu+ mu-", "eta'", "pi+,pi-,mu+,mu-", 2.9E-5);     // upper limit
    sdata->AddDecay("eta' --> pi+ pi- pi0",     "eta'", "pi+,pi-,pi0",     3.8E-4);
    sdata->AddDecay("eta' --> pi+ pi- e+ e-",   "eta'", "pi+,pi-,e+,e-",   2.4E-3);

    std::cout << "=============================================" << endl;
    std::cout << "eta' after change:" << endl;
    listParticle("eta'");

    std::cout << "=============================================" << endl;
    std::cout << "rho0 before change:" << endl;
    listParticle("rho0");
    // additional rho decays ( PDG 2016)
    sdata->AddDecay("rho0 --> pi+ pi- g",           "rho0", "pi+,pi-,g",            9.9E-3);
    sdata->AddDecay("rho0 --> pi0 g",               "rho0", "pi0,g",                6.0E-4);
    sdata->AddDecay("rho0 --> eta g",               "rho0", "eta,g",                3.0E-4);
//    sdata->AddDecay("rho0 --> eta g",               "rho0", "eta,g",                3.0E-4);
    sdata->AddDecay("rho0 --> pi0 pi0 g",           "rho0", "pi0,pi0,g",            4.5E-5);
    sdata->AddDecay("rho0 --> pi+ pi- pi0 pi0 g",   "rho0", "pi+,pi-,pi0,pi0,g",    1.6E-5);
    sdata->AddDecay("rho0 --> pi+ pi- pi+ pi- g",   "rho0", "pi+,pi-,pi+,pi-,g",    1.8E-5);
    std::cout << "=============================================" << endl;
    std::cout << "rho0 after change:" << endl;
    listParticle("rho0");
}

}
}
}
