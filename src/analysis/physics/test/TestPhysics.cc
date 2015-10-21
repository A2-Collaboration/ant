#include "TestPhysics.h"

#include "data/Event.h"
#include "data/Particle.h"
#include "base/ParticleType.h"
#include "utils/combinatorics.h"
#include "TH1D.h"
#include <memory>
#include <iostream>
#include "plot/HistogramFactories.h"
#include "TCanvas.h"

#include "plot/root_draw.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace ant::analysis::data;

ParticleCombinatoricsTest::ParticleCombinatoricsTest(const std::string& name, PhysOptPtr opts):
    Physics(name, opts)
{
    const BinSettings im_binning(100,0,250);
    const BinSettings energy_binning(100,0,250);
    const BinSettings npart_binning(10,0,10);

    ggim     = HistFac.makeHist<double>("2 #gamma IM", "M_{#gamma #gamma} [MeV]","#", im_binning);
    gggim    = HistFac.makeHist<double>("3 #gamma im","M_{#gamma #gamma #gamma} [MeV]","#", im_binning);
    nphotons = HistFac.makeHist<int>("Number of photons", "N", "", npart_binning);
    nprotons = HistFac.makeHist<int>("Number of protons","N","",npart_binning);

    // Build a map of ParticleType -> Histogram, and fill it
    for( auto& type : ParticleTypeDatabase() ) {
        EHists[&type] = HistFac.KinEnergyPlot(type.PrintName()+" Energy");
    }

}


void ParticleCombinatoricsTest::ProcessEvent(const Event &event)

{

    const ParticleList& photons = event.Reconstructed().Particles().Get(ParticleTypeDatabase::Photon);
    const ParticleList& protons = event.Reconstructed().Particles().Get(ParticleTypeDatabase::Proton);
    const ParticleList& all = event.Reconstructed().Particles().GetAll();

    for( auto& particle : all ) {

        // fill the histogram corresponding to the partice type of the current particle
        auto entry = EHists.find(&particle->Type());
        if( entry != EHists.end()) {
            entry->second.Fill(particle);
        }

    }

    nphotons.Fill(photons.size());
    nprotons.Fill(protons.size());

    auto combinations2 = utils::makeCombination(photons,2);
    do {
        TLorentzVector v;
        for( auto& i: combinations2 ) {
            v += *i;
        }

        ggim.Fill(v.M());

    } while(combinations2.next());

    auto combinations3 = utils::makeCombination(photons,3);
    do {
        TLorentzVector v;
        for( auto& i: combinations3 ) {
            v += *i;
        }

        gggim.Fill(v.M());

    } while(combinations3.next());
}

void ParticleCombinatoricsTest::Finish()
{

}

void ParticleCombinatoricsTest::ShowResult()
{
    canvas cc("ParticleCombinatoricsTest");
    cc << ggim << gggim << nphotons << nprotons << endc;

    canvas cc2("ParticleCombinatoricsTest2");
    for( auto& e : EHists) {
        cc2 << e.second;
    }
    cc2 << endc;

}

AUTO_REGISTER_PHYSICS(ParticleCombinatoricsTest)
