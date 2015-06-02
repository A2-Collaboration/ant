#include "omega_bottomup.h"

#include "plot/root_draw.h"
#include "plot/HistogramFactories.h"

#include <algorithm>
#include <iostream>

using namespace ant;
using namespace std;

void ant::analysis::OmegaBottomUp::findPossibleDecays(const ParticlePtr gamma, const ParticlePtr meson_gamma_1, const ParticlePtr meson_gamma_2, decaylist_t& decays)
{

    TLorentzVector meson_cand = TLorentzVector(*meson_gamma_1) + TLorentzVector(*meson_gamma_2);
    const ParticleTypeDatabase::Type* t = nullptr;

    if( eta_mass_cut.Contains(meson_cand.M())) {
        t = &ParticleTypeDatabase::Eta;
    } else if( pi0_mass_cut.Contains(meson_cand.M())) {
        t = &ParticleTypeDatabase::Pi0;
    }

    if(t) {
        TLorentzVector omega_cand = meson_cand + TLorentzVector(*gamma);

        if(omega_mass_cut.Contains(omega_cand.M())) {
            decays.emplace_back( omega_decay(
                ParticlePtr(new Particle(ParticleTypeDatabase::Omega, omega_cand)),
                ParticlePtr(new Particle(*t, meson_cand))));
        }
    }
}

ant::analysis::OmegaBottomUp::OmegaBottomUp(const std::string& name):
    Physics(name),
    eta_mass_cut(IntervalD::CenterWidth(ParticleTypeDatabase::Eta.Mass(), 50.0)),
    pi0_mass_cut(IntervalD::CenterWidth(ParticleTypeDatabase::Pi0.Mass(), 20.0)),
    omega_mass_cut(IntervalD::CenterWidth(ParticleTypeDatabase::Omega.Mass(), 80.0))
{
    omega_eta_found = HistFac.makeHist<int>(
                "#omega #rightarrow #eta #gamma #rightarrow (#gamma #gamma) #gamma per event",
                "number of decays / event",
                "",
                BinSettings(10),
                "omega_eta_per_event");

    omega_pi0_found = HistFac.makeHist<int>(
                "#omega #rightarrow #pi^{0} #gamma #rightarrow (#gamma #gamma) #gamma per event",
                "number of decays / event",
                "",
                BinSettings(10),
                "omega_pi0_per_event");

    omega_IM = HistFac.InvariantMass("#omega IM");
    eta_IM   = HistFac.InvariantMass("eta IM");
    pi0_IM   = HistFac.InvariantMass(ParticleTypeDatabase::Pi0.PrintName()+" IM");

}

void ant::analysis::OmegaBottomUp::ProcessEvent(const ant::Event &event)
{

    decaylist_t decays;

    ParticleList photons = event.Reconstructed().Particles().Get(ParticleTypeDatabase::Photon);

    sort(photons.begin(),photons.end(),[] (const ParticlePtr& a, const ParticlePtr& b) { return a->E() > b->E();});

    if(photons.size() < 3) return;

    findPossibleDecays(photons.at(0), photons.at(1), photons.at(2), decays);
    findPossibleDecays(photons.at(1), photons.at(0), photons.at(2), decays);
    findPossibleDecays(photons.at(2), photons.at(0), photons.at(1), decays);



    int n_omega_eta =0;
    int n_omega_pi0 =0;
    for(auto& d : decays) {
        omega_IM.Fill(d.omega);
        if(d.meson2->Type() == ParticleTypeDatabase::Eta) {
            eta_IM.Fill(d.meson2);
            n_omega_eta++;
        } else if(d.meson2->Type() == ParticleTypeDatabase::Pi0) {
            pi0_IM.Fill(d.meson2);
            n_omega_pi0++;
        }
    }
    omega_eta_found.Fill(n_omega_eta);
    omega_pi0_found.Fill(n_omega_pi0);
}

void ant::analysis::OmegaBottomUp::Finish()
{

}

void ant::analysis::OmegaBottomUp::ShowResult()
{
    canvas("Omega Buttom Up") << omega_eta_found << omega_pi0_found << omega_IM << eta_IM << pi0_IM << endc;
}
