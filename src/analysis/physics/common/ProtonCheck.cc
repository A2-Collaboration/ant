#include "physics/common/ProtonCheck.h"
#include "plot/root_draw.h"
#include "data/Particle.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include <cmath>
#include <iostream>

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace ant::analysis::data;

ProtonCheck::ProtonCheck(const std::string& name,PhysOptPtr opts):
    Physics(name,opts)
{
    const BinSettings e(1000);
    const BinSettings t(300,-15,15);
    const BinSettings dE(80,0,8);
    const BinSettings theta_bins(180,0,180);
    tof = HistFac.makeTH2D("Proton TOF TAPS","t [ns]","E [MeV]", t,e,"tof");
    tof_trueE = HistFac.makeTH2D("Proton TOF TAPS (MCTrue Energy)","t [ns]","E_{true} [MeV]", t,e,"tof_true");
    dEE = HistFac.makeTH2D("Proton dEE TAPS","E [MeV]","dE [MeV]", e,dE,"dEE");
    cand_mult = HistFac.makeTH1D("Candidates / Event","# Candiadates/Event","#",BinSettings(20),"mult");

    theta =  HistFac.makeTH1D("Theta","#Theta","#",BinSettings(180,0,180),"theta");

    theta_corr = HistFac.makeTH2D("Theta Corrleation","true","rec",theta_bins,theta_bins,"theta_corr");
}

void ProtonCheck::ProcessEvent(const Event &event)
{
    if(event.MCTrue().Particles().GetAll().size() == 1) {
        const auto& protons = event.MCTrue().Particles().Get(ParticleTypeDatabase::Proton);

        if(protons.size()==1) {
            const auto& mctrue = protons.at(0);

            if(mctrue->Theta() < 20.0*TMath::DegToRad()) {

                for(const auto& cand : event.Reconstructed().Candidates()) {

                    if(cand->Detector() == Detector_t::Any_t::TAPS) {
                        tof->Fill(cand->Time(), cand->ClusterEnergy());
                        tof_trueE->Fill(cand->Time(), mctrue->Ek());
                        dEE->Fill(cand->ClusterEnergy(), cand->VetoEnergy());
                        theta->Fill(cand->Theta()*TMath::RadToDeg());
                        theta_corr->Fill(mctrue->Theta()*TMath::RadToDeg(), cand->Theta()*TMath::RadToDeg());
                    }
                }

                cand_mult->Fill(event.Reconstructed().Candidates().size());
            }

        }

    }
}


void ProtonCheck::Finish()
{

}


void ProtonCheck::ShowResult()
{
    canvas("ProtonCheck")
            << cand_mult
            << drawoption("colz")
            << tof << tof_trueE << dEE << theta_corr << endc;
}


AUTO_REGISTER_PHYSICS(ProtonCheck)
