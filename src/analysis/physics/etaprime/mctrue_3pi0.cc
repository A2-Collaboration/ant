#include "mctrue_3pi0.h"
#include "base/std_ext/math.h"
#include "base/std_ext/string.h"
#include "utils/ParticleTools.h"

#include "TTree.h"
#include "TCanvas.h"

using namespace std;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::utils;
using namespace ant::analysis::physics;



std::vector<double> McTrue3Pi0::GetAllPhotonAngles(const TParticleList& photons) const
{
    vector<double> angles((photons.size() * (photons.size() - 1)) / 2);
    unsigned index = 0;
    for ( unsigned i = 0 ; i < photons.size() ; ++i )
        for ( unsigned j = i + 1 ; j < photons.size() ; ++j)
        {
            vec3 ph_i(photons.at(i)->p);
            vec3 ph_j(photons.at(j)->p);

            angles.at(index) = ph_i.Angle(ph_j);
            index++;
        }
    return angles;
}

McTrue3Pi0::McTrue3Pi0(const std::string& name, OptionsPtr opts) :
    Physics(name, opts),
    pi0s(3),
    popens(3),
    gammas(6)
{
    mcTrue = HistFac.makeTTree("mcTrue");

    proton.SetBranches(mcTrue,"proton");
    for (int i = 0 ; i < 3 ; ++i)
        pi0s.at(i).SetBranches(mcTrue,formatter() << "pi0_" << i );
    for (int i = 0 ; i < 6 ; ++i)
        gammas.at(i).SetBranches(mcTrue,formatter() << "gamma_" << i );
    mcTrue->Branch("p0open",addressof(popens.at(0)));
    mcTrue->Branch("p1open",addressof(popens.at(1)));
    mcTrue->Branch("p2open",addressof(popens.at(2)));

    hAngle = HistFac.makeTH1D("Photon Angles","#alpha [#circ]","#",BinSettings(180));
}


void McTrue3Pi0::ProcessEvent(const TEvent& event, manager_t&)
{
    const auto& mcdata = event.MCTrue();

    /// \todo this could actually be implemented with ParticleTree_t comparison

    const auto& pions  = utils::ParticleTools::FindParticles(ParticleTypeDatabase::Pi0, mcdata.ParticleTree);
    auto mctrue_particles = utils::ParticleTypeList::Make(event.MCTrue().ParticleTree);
    const auto& protons = mctrue_particles.Get(ParticleTypeDatabase::Proton);
    const auto& photons = mctrue_particles.Get(ParticleTypeDatabase::Photon);

    if (pions.size() == 3 && protons.size() == 1 && photons.size() == 6)
    {
        proton = ParticleVars(*(protons.at(0)));
        for (int i = 0 ; i < 3 ; ++i)
        {
            pi0s.at(i) = ParticleVars(*(pions.at(i)));
            popens.at(i) = std_ext::radian_to_degree(photons.at(2*i)->p.Angle(photons.at(2*i+1)->p));
        }
        for (int i = 0 ; i < 6 ; ++i)
            gammas.at(i) = ParticleVars(*photons.at(i));
    }

    for ( const auto& angle: GetAllPhotonAngles(photons))
        hAngle->Fill(std_ext::radian_to_degree(angle));

    mcTrue->Fill();
}

void McTrue3Pi0::ShowResult()
{
    canvas("Angles") << hAngle << endc;
}




AUTO_REGISTER_PHYSICS(McTrue3Pi0)
