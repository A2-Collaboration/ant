#include "PhysicsStats.h"
#include "data/Particle.h"
#include "base/ParticleType.h"
#include "plot/root_draw.h"
#include <sstream>

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace ant::analysis::data;


TPC_PhysicsStats::TPC_PhysicsStats(PhysOptPtr& opts): Physics("TPC_PhysicsStats", opts)
{

}


std::string ProdProcess(const ParticlePtr& beamtarget) {

    stringstream s;

//    if(beamtarget->Type() == ParticleTypeDatabase::BeamProton || beamtarget->Type() == ParticleTypeDatabase::BeamNeutron) {
//        s << beamtarget->Type().PrintName() << " #rightarrow ";
//        for(const auto& d : beamtarget->Daughters()) {
//            s << " " << d->Type().PrintName();
//        }
//    }

    return s.str();
}

bool containsCharged(const ParticleList& particles) {
    for(const auto& p : particles) {
        if(p->Type() == ParticleTypeDatabase::PiCharged || p->Type() == ParticleTypeDatabase::eCharged)
            return true;
    }
    return false;
}

void TPC_PhysicsStats::ProcessEvent(const Event& event)
{
    const auto& intermediates = event.MCTrue().Intermediates().GetAll();
    if(!intermediates.empty()) {
        const auto process = ProdProcess(intermediates.front());
        TH1D* h = nullptr;

        const auto entry = decay_angles.find(process);

        if(entry == decay_angles.end()) {
            h = HistFac.makeTH1D(process,"#Theta [#circ]","",BinSettings(360,0,180),process+"_ch_angle");
            decay_angles.insert({process, h});
        } else {
            h = entry->second;
        }

        const auto& particles = event.MCTrue().Particles().GetAll();
        if(containsCharged(particles)) {

            for(const ParticlePtr& p : particles) {
                if(p->Type().Charged()) {
                    h->Fill(p->Theta()*TMath::RadToDeg());
                }
            }
        }

    }
}

void TPC_PhysicsStats::Finish()
{

}

void TPC_PhysicsStats::ShowResult()
{
    canvas c("TPC_PhysicsStats");
    for(const auto& e : decay_angles) {
       c << e.second;
    }

    c << endc;
}

AUTO_REGISTER_PHYSICS(TPC_PhysicsStats, "TPC_PhysicsStats")
