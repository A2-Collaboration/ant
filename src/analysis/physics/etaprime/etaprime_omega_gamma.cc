#include "etaprime_omega_gamma.h"
#include "plot/root_draw.h"

using namespace std;
using namespace ant::analysis;
using namespace ant::analysis::physics;


EtapOmegaG::EtapOmegaG(PhysOptPtr opts) : Physics("EtapOmegaG", opts)
{
    gggg = HistFac.makeTH1D("gggg","4#gamma IM","events",bins_im,"gggg");
}




void ant::analysis::physics::EtapOmegaG::ProcessEvent(const data::Event& event)
{
    const auto& data = event.MCTrue();

    const auto nParticles = data.Particles().GetAll().size();
    if(nParticles != 4 && nParticles != 5)
        return;

    const auto& photons = data.Particles().Get(ParticleTypeDatabase::Photon);
    const auto& protons = data.Particles().Get(ParticleTypeDatabase::Proton);

    const auto nPhotons = photons.size();
    const auto nProtons = protons.size();

    if(nPhotons != 4 && nProtons+nPhotons != nParticles)
        return;

    TLorentzVector sum(0,0,0,0);
    for(const auto& photon : photons) {
        sum += *photon;
    }
    gggg->Fill(sum.M());
}

void ant::analysis::physics::EtapOmegaG::ShowResult()
{
    canvas(GetName()) << gggg << endc;
}


AUTO_REGISTER_PHYSICS(EtapOmegaG, "EtapOmegaG")
