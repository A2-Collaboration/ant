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

void EtapOmegaG::ProcessEvent(const data::Event& event)
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
    auto fill_combinations = [] (TH1* h, unsigned multiplicity, const data::ParticleList& particles) {
        for( auto comb = utils::makeCombination(particles,multiplicity); !comb.Done(); ++comb) {
             TLorentzVector sum(0,0,0,0);
             for(const auto& p : comb) {
                 sum += *p;
             }
             h->Fill(sum.M());
        }
    };

    fill_combinations(gg,   2, photons);
    fill_combinations(ggg,  3, photons);
    fill_combinations(gggg, 4, photons);

}

void EtapOmegaG::ShowResult()
{
    canvas(GetName()) << gg << ggg << gggg << endc;
}


AUTO_REGISTER_PHYSICS(EtapOmegaG, "EtapOmegaG")
