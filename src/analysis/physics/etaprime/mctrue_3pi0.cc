#include "mctrue_3pi0.h"
#include "plot/root_draw.h"
#include "base/std_ext/math.h"
#include "base/std_ext/string.h"
#include "data/Particle.h"
#include "utils/particle_tools.h"

#include "TTree.h"
#include "TCanvas.h"

using namespace std;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::utils;
using namespace ant::analysis::data;
using namespace ant::analysis::physics;



std::vector<double> McTrue3Pi0::GetAllPhotonAngles(const ParticleList& photons) const
{
    vector<double> angles((photons.size() * (photons.size() - 1)) / 2);
    for ( unsigned i = 0 ; i < photons.size() ; ++i )
        for ( unsigned j = i + 1 ; j < photons.size() ; ++j)
        {
            TVector3 ph_i(photons.at(i)->Vect());
            TVector3 ph_j(photons.at(j)->Vect());

            angles.push_back(ph_i.Angle(ph_j));
        }
    return angles;
}

McTrue3Pi0::McTrue3Pi0(const std::string& name, PhysOptPtr opts) :
    Physics(name, opts),
    pi0s(3),
    gammas(6)
{
    mcTrue = HistFac.makeTTree("mcTrue");

    proton.SetBranches(mcTrue,"proton");
    for (int i = 0 ; i < 3 ; ++i)
        pi0s.at(i).SetBranches(mcTrue,formatter() << "pi0_" << i );
    for (int i = 0 ; i < 6 ; ++i)
        gammas.at(i).SetBranches(mcTrue,formatter() << "gamma_" << i );

    hAngle = HistFac.makeTH1D("Photon Angles","#alpha [#circ]","#",BinSettings(100,0,30));
}


void McTrue3Pi0::ProcessEvent(const data::Event& event)
{
    const auto& mcdata = event.MCTrue();

    const auto& pions  = mcdata.Intermediates().Get(ParticleTypeDatabase::Pi0);
    const auto& protons = mcdata.Particles().Get(ParticleTypeDatabase::Proton);
    const auto& photons = mcdata.Particles().Get(ParticleTypeDatabase::Photon);

    if (pions.size() == 3 && protons.size() == 1 && photons.size() == 6)
    {
        proton = ParticleVars(*(protons.at(0)));
        for (int i = 0 ; i < 3 ; ++i)
            pi0s.at(i) = ParticleVars(*(pions.at(i)));
        for (int i = 0 ; i < 6 ; ++i)
            gammas.at(i) = ParticleVars(*photons.at(i));
    }

    for ( const auto& angle: GetAllPhotonAngles(photons))
        hAngle->Fill(angle * TMath::RadToDeg());

    mcTrue->Fill();
}

void McTrue3Pi0::ShowResult()
{
    canvas("Angles") << hAngle << endc;
}




AUTO_REGISTER_PHYSICS(McTrue3Pi0)
