#include "MCSmearing.h"

#include "base/std_ext/math.h"

#include "TDirectory.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

MCSmearing::MCSmearing(PhysOptPtr opts) : Physics("MCSmearing", opts)
{
    BinSettings bins_energies(50, 0, 1000);
    BinSettings bins_theta(50, 0, 180);

    angles = HistFac.makeTH3D("Rec-True Angles","True Energy / MeV","True #theta / #circ","Difference Angle / #circ",
                              bins_energies,
                              bins_theta,
                              BinSettings(100, 0, 15),
                              "angles"
                              );
    energies = HistFac.makeTH3D("(Rec-True)/True Energies","True Energy / MeV","True #theta / #circ","Energies / #circ",
                              bins_energies,
                              bins_theta,
                              BinSettings(100, -0.3, 0.05),
                              "energies"
                              );
}

MCSmearing::~MCSmearing()
{

}

void MCSmearing::ProcessEvent(const data::Event& event)
{

    const auto& true_photons = event.MCTrue().Particles().Get(ParticleTypeDatabase::Photon);
    const auto& reco_photons = event.Reconstructed().Particles().Get(ParticleTypeDatabase::Photon);

    if(true_photons.size() != 1 || reco_photons.size() != 1)
        return;

    const data::ParticlePtr& true_g = true_photons.front();
    const data::ParticlePtr& reco_g = reco_photons.front();

    const auto& true_theta = std_ext::radian_to_degree(true_g->Theta());


    angles->Fill(true_g->Ek(), true_theta, std_ext::radian_to_degree(reco_g->Angle(true_g->Vect())));
    energies->Fill(true_g->Ek(), true_theta, (reco_g->Ek()-true_g->Ek())/true_g->Ek());


}

void MCSmearing::Finish()
{
    energies->FitSlicesZ();
    angles->FitSlicesZ();
}

void MCSmearing::ShowResult()
{
    canvas("MCSmearing") << drawoption("colz")
                         << (TH2D*)gDirectory->Get("energies_1")
                         << (TH2D*)gDirectory->Get("energies_2")
                         << (TH2D*)gDirectory->Get("angles_1")
                         << (TH2D*)gDirectory->Get("angles_2")
                         << endc;
}


AUTO_REGISTER_PHYSICS(MCSmearing, "MCSmearing")
