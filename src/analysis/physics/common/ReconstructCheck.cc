#include "physics/common/ReconstructCheck.h"
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

ReconstructCheck::ReconstructCheck(PhysOptPtr opts):
    Physics("ReconstructCheck",opts)
{
    const BinSettings e(max(1000.0, atof(GetOption("Emax").c_str())));
    EnergyRec_cb = HistFac.makeTH2D("Energry Reconstruction CB","E_{true} [MeV]","E_{rec} [MeV]", e, e, "Energy_rec_cb");
    EnergyRec_taps = HistFac.makeTH2D("Energry Reconstruction TAPS","E_{true} [MeV]","E_{rec} [MeV]", e, e, "Energy_rec_taps");
}

void ReconstructCheck::ProcessEvent(const Event &event)
{
    if(event.MCTrue().Particles().GetAll().size() == 1) {
        const auto& p = event.MCTrue().Particles().GetAll().at(0);

        if( event.Reconstructed().Candidates().size() == 1)
        for(auto& cand : event.Reconstructed().Candidates()) {

            if(cand->Detector() & Detector_t::Any_t::CB) {
                EnergyRec_cb->Fill(p->Ek(), cand->ClusterEnergy());
            } else if( cand->Detector() & Detector_t::Any_t::TAPS) {
                EnergyRec_taps->Fill(p->Ek(), cand->ClusterEnergy());
            }

        }

    }
}


void ReconstructCheck::Finish()
{

}


void ReconstructCheck::ShowResult()
{
    canvas("ReconstructCheck")
            << drawoption("colz")
            << EnergyRec_cb << EnergyRec_taps << endc;
}


AUTO_REGISTER_PHYSICS(ReconstructCheck, "ReconstructCheck")
