#include "etaprime_omega_gamma.h"
#include "plot/root_draw.h"
#include "utils/combinatorics.h"

using namespace std;
using namespace ant::analysis;
using namespace ant::analysis::physics;


EtapOmegaG::EtapOmegaG(PhysOptPtr opts) : Physics("EtapOmegaG", opts)
{
    gggg = HistFac.makeTH1D("gggg","4#gamma IM / MeV","events",bins_im,"gggg");
    ggg = HistFac.makeTH1D("ggg","3#gamma IM / MeV","events",bins_im,"ggg");
    gg = HistFac.makeTH1D("gg","2#gamma IM / MeV","events",bins_im,"gg");
}





void ant::analysis::physics::EtapOmegaG::ProcessEvent(const data::Event& event)
{
    const auto& data = event.Reconstructed();

    const auto nParticles = data.Particles().GetAll().size();
    if(nParticles != 4 && nParticles != 5)
        return;

    const auto& photons = data.Particles().Get(ParticleTypeDatabase::Photon);
    const auto& protons = data.Particles().Get(ParticleTypeDatabase::Proton);

    const auto nPhotons = photons.size();
    const auto nProtons = protons.size();

    if(nPhotons != 4 && nProtons+nPhotons != nParticles)
        return;



    // gamma combinatorics


    for( auto gcomb = utils::makeCombination(photons,2); !gcomb.Done(); ++gcomb) {
         TLorentzVector sum2(0,0,0,0);
         for(const auto& photon : gcomb) {
             sum2 += *photon;
         }
         gg->Fill(sum2.M());
    }

    for( auto gcomb = utils::makeCombination(photons,3); !gcomb.Done(); ++gcomb) {
         TLorentzVector sum3(0,0,0,0);
         for(const auto& photon : gcomb) {
             sum3 += *photon;
         }
         ggg->Fill(sum3.M());
    }

    TLorentzVector sum4(0,0,0,0);
    for(const auto& photon : photons) {
        sum4 += *photon;
    }
    gggg->Fill(sum4.M());

}

void ant::analysis::physics::EtapOmegaG::ShowResult()
{
    canvas(GetName()) << gg << ggg << gggg << endc;
}


AUTO_REGISTER_PHYSICS(EtapOmegaG, "EtapOmegaG")
