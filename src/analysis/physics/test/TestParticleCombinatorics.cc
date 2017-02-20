#include "TestParticleCombinatorics.h"

#include "base/ParticleType.h"
#include "utils/combinatorics.h"
#include "TH1D.h"
#include <memory>
#include <iostream>
#include "plot/HistogramFactories.h"
#include "TCanvas.h"

#include "plot/root_draw.h"

#include "utils/particle_tools.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

TestParticleCombinatorics::TestParticleCombinatorics(const std::string& name, OptionsPtr opts):
    Physics(name, opts)
{
    const BinSettings im_binning(100,0,250);
    const BinSettings npart_binning(10,0,10);
    const BinSettings energy_binning(100,0,250);

    ggim     = HistFac.makeTH1D("2 #gamma IM",       "M_{#gamma #gamma} [MeV]",       "#", im_binning);
    gggim    = HistFac.makeTH1D("3 #gamma im",       "M_{#gamma #gamma #gamma} [MeV]","#", im_binning);
    nphotons = HistFac.makeTH1D("Number of photons", "N",                             "",  npart_binning);
    nprotons = HistFac.makeTH1D("Number of protons", "N",                             "",  npart_binning);

    // Build a map of ParticleType -> Histogram, and fill it
    for( auto& type : ParticleTypeDatabase() ) {
        EHists[&type] = HistFac.makeTH1D(type.PrintName()+" Energy", "E [MeV]", "", energy_binning);
    }

}


void TestParticleCombinatorics::ProcessEvent(const TEvent& event, manager_t&)
{
    auto recon_particles = utils::ParticleTypeList::Make(event.Reconstructed().Candidates);

    const TParticleList& photons = recon_particles.Get(ParticleTypeDatabase::Photon);
    const TParticleList& protons = recon_particles.Get(ParticleTypeDatabase::Proton);
    const TParticleList& all = recon_particles.GetAll();

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
        LorentzVec v;
        for( auto& i: combinations2 ) {
            v += *i;
        }

        ggim->Fill(v.M());

    } while(combinations2.next());

    auto combinations3 = utils::makeCombination(photons,3);
    do {
        LorentzVec v;
        for( auto& i: combinations3 ) {
            v += *i;
        }

        gggim->Fill(v.M());

    } while(combinations3.next());
}

void TestParticleCombinatorics::Finish()
{

}

void TestParticleCombinatorics::ShowResult()
{
    canvas cc("ParticleCombinatoricsTest");
    cc << ggim << gggim << nphotons << nprotons << endc;

    canvas cc2("ParticleCombinatoricsTest2");
    for( auto& e : EHists) {
        cc2 << e.second;
    }
    cc2 << endc;

}

AUTO_REGISTER_PHYSICS(TestParticleCombinatorics)
