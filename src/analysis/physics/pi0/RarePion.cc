#include "RarePion.h"

#include "base/ParticleType.h"
#include "utils/combinatorics.h"
#include "TH1D.h"
#include <memory>
#include <iostream>
#include "plot/HistogramFactories.h"
#include "TCanvas.h"

#include "plot/root_draw.h"

#include "utils/particle_tools.h"

/**
 * this is mainly a copy of the class TestPhysics and to be changed that it fits the analysis
 * of rare pion decays.
 */

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

RarePion::RarePion(const std::string& name, OptionsPtr opts):
    Physics(name, opts)
{
    const BinSettings im_binning(100,0,250);
    const BinSettings npart_binning(10,0,10);
    const BinSettings energy_binning(100,0,250);

    ggim     = HistFac.makeTH1D("2 #gamma IM",       "M_{#gamma #gamma} [MeV]",       "#", im_binning);
    gggim    = HistFac.makeTH1D("3 #gamma im",       "M_{#gamma #gamma #gamma} [MeV]","#", im_binning);
    ggggim   = HistFac.makeTH1D("4 #gamma im",       "M_{#gamma #gamma #gamma #gamma} [MeV]","#", im_binning);
    nphotons = HistFac.makeTH1D("Number of photons", "N",                             "",  npart_binning);
    nprotons = HistFac.makeTH1D("Number of protons", "N",                             "",  npart_binning);

    // Build a map of ParticleType -> Histogram, and fill it
    for( auto& type : ParticleTypeDatabase() ) {
        EHists[&type] = HistFac.makeTH1D(type.PrintName()+" Energy", "E [MeV]", "", energy_binning);
    }

}


void RarePion::ProcessEvent(const TEvent& event, manager_t&)

{
    const TParticleList& photons = event.Reconstructed().Particles.Get(ParticleTypeDatabase::Photon);
    const TParticleList& protons = event.Reconstructed().Particles.Get(ParticleTypeDatabase::Proton);
    const TParticleList& all = event.Reconstructed().Particles.GetAll();

    for( auto& particle : all ) {

        // fill the histogram corresponding to the partice type of the current particle
        auto entry = EHists.find(&particle->Type());
        if( entry != EHists.end()) {
            entry->second->Fill(particle->Ek());
        }

    }

    nphotons->Fill(photons.size());
    nprotons->Fill(protons.size());

    auto combinations2 = utils::makeCombination(photons,2);
    do {
        TLorentzVector v;
        for( auto& i: combinations2 ) {
            v += *i;
        }

        ggim->Fill(v.M());

    } while(combinations2.next());

    auto combinations3 = utils::makeCombination(photons,3);
    do {
        TLorentzVector v;
        for( auto& i: combinations3 ) {
            v += *i;
        }

        gggim->Fill(v.M());

    } while(combinations3.next());

    auto combinations4 = utils::makeCombination(photons,4);
    do {
        TLorentzVector v;
        for( auto& i: combinations4 ) {
            v += *i;
        }

        ggggim->Fill(v.M());

    } while(combinations4.next());
}

void RarePion::Finish()
{

}

void RarePion::ShowResult()
{
    canvas cc("RarePion");
    cc <<
          ggim <<
          gggim <<
          ggggim <<
          nphotons <<
          //nprotons <<
          endc;

//    canvas cc2("RarePion2");
//    for( auto& e : EHists) {
//        cc2 << e.second;
//    }
//    cc2 << endc;

}

AUTO_REGISTER_PHYSICS(RarePion)
