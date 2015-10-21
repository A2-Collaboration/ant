#include "MCSmearing.h"

#include "base/std_ext/math.h"

#include "utils/particle_tools.h"

#include "TDirectory.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

MCSmearing::MCSmearing(const std::string& name, PhysOptPtr opts) : Physics(name, opts)
{
    BinSettings bins_energies(50, 0, 1000);
    BinSettings bins_theta(50, 0, 180);
    BinSettings bins_im(2000, 0, 1000);

    angles = HistFac.makeTH3D("Rec-True Angles","True Energy / MeV","True #theta / #circ","Difference Angle / #circ",
                              bins_energies,
                              bins_theta,
                              BinSettings(100, 0, 15),
                              "angles"
                              );
    energies = HistFac.makeTH3D("(Rec-True)/True Energies","True Energy / MeV","True #theta / #circ","Energies / MeV",
                              bins_energies,
                              bins_theta,
                              BinSettings(100, -0.3, 0.05),
                              "energies"
                              );
    IM = HistFac.makeTH1D("Invariant mass","IM / MeV","#",
                          bins_im,
                          "IM"
                          );
    IM_2g = HistFac.makeTH1D("Invariant mass 2#gamma","IM / MeV","#",
                          bins_im,
                          "IM_2g"
                          );
}

MCSmearing::~MCSmearing()
{

}

void MCSmearing::ProcessEvent(const data::Event& event)
{

    const auto& true_photons = event.MCTrue().Particles().Get(ParticleTypeDatabase::Photon);
    const auto& reco_photons = event.Reconstructed().Particles().Get(ParticleTypeDatabase::Photon);

    utils::ParticleTools::FillIMCombinations(IM, reco_photons.size(), reco_photons);

    if(true_photons.size() == 1 && reco_photons.size() == 1) {

        const data::ParticlePtr& true_g = true_photons.front();
        const data::ParticlePtr& reco_g = reco_photons.front();

        const auto& true_theta = std_ext::radian_to_degree(true_g->Theta());

        angles->Fill(true_g->Ek(), true_theta, std_ext::radian_to_degree(reco_g->Angle(true_g->Vect())));
        energies->Fill(true_g->Ek(), true_theta, (reco_g->Ek()-true_g->Ek())/true_g->Ek());
    }

    if(true_photons.size() == 2 && reco_photons.size() == 2) {
        const TLorentzVector sum(*reco_photons.front() + *reco_photons.back());
        IM_2g->Fill(sum.M());
    }

}

void MCSmearing::Finish()
{
    if(energies->GetEntries()>0)
        energies->FitSlicesZ();
    if(angles->GetEntries()>0)
        angles->FitSlicesZ();
}

void MCSmearing::ShowResult()
{
    canvas c(GetName());

    c << IM << IM_2g;

    if(energies->GetEntries()>0)
        c << drawoption("colz")
          << (TH2D*)gDirectory->Get("energies_1")
          << (TH2D*)gDirectory->Get("energies_2");
    if(angles->GetEntries()>0)
        c << (TH2D*)gDirectory->Get("angles_1")
          << (TH2D*)gDirectory->Get("angles_2");

    c << endc;
}


AUTO_REGISTER_PHYSICS(MCSmearing)
