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



McTrue3Pi0::McTrue3Pi0(const std::string& name, PhysOptPtr opts) :
    Physics(name, opts),
    pi0s(3)
{
    mcTrue = HistFac.makeTTree("mcTrue");

    proton.SetBranches(mcTrue,"proton");
    for (int i = 0 ; i < 3 ; ++i)
        pi0s.at(i).SetBranches(mcTrue,formatter() << "pi0_" << i );

}


void McTrue3Pi0::ProcessEvent(const data::Event& event)
{
    const auto& mcdata = event.MCTrue();

    const auto& pions  = mcdata.Intermediates().Get(ParticleTypeDatabase::Pi0);
    const auto& protons = mcdata.Particles().Get(ParticleTypeDatabase::Proton);

    if (pions.size() == 3 && protons.size() == 1)
    {
        proton = ParticleVars(*(protons.at(0)));
        for (int i = 0 ; i < 3 ; ++i)
            pi0s.at(i) = ParticleVars(*(pions.at(i)));
    }

    mcTrue->Fill();
//    h6photonEvents->Fill(utils::ParticleTools::GetDecayString(mcdata.ParticleTree()).c_str(),1);
}

void McTrue3Pi0::ShowResult()
{
}




AUTO_REGISTER_PHYSICS(McTrue3Pi0)
