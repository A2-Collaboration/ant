#include "MCOverview.h"
#include "Particle.h"
#include "Track.h"
#include "utils/combinatorics.h"
#include "plot/root_draw.h"
#include "Event.h"
#include "TMath.h"
#include <tuple>
#include <iostream>

using namespace std;
using namespace ant;

const ant::ParticleTypeDatabase::Type* GetParticleType(const ParticlePtr& p )
{
    return &p->Type();
}

ant::analysis::MCOverview::MCOverview(const mev_t energy_scale)
{
    HistogramFactory::SetName("MCOverview");

    const BinSettings energy_bins(100,0,energy_scale);
    const BinSettings theta_bins(100,0,180);
    const BinSettings theta_proton_bins(100,0,60.0);
    const BinSettings particle_type_bins(10,0,10);

    mc_particle_stats.AddHist1D(
                [] (const ParticlePtr& p) { return p->Ek();},
                HistogramFactory::Make1D("MC Energy","E [MeV]","", energy_bins)
    );

    mc_particle_stats.AddHist1D(
                [] (const ParticlePtr& p) { return p->Type().PrintName().c_str();},
                HistogramFactory::Make1D(
                    "Generated Particles",
                    "Particle Type",
                    "#",
                    particle_type_bins,
                    "ParticleTypes"
                    )
    );

    mc_particle_stats.AddHist2D(
                [] (const ParticlePtr& p) { return make_tuple(p->Ek(), p->Theta()*TMath::RadToDeg());},
                HistogramFactory::Make2D(
            "MC Energy",
            "E [MeV]",
            "#theta [#circ]",
            energy_bins,
            theta_bins)
    );

    auto ptype = mc_particle_stats.AddBranchNode<const ParticleTypeDatabase::Type*>(GetParticleType);

    for( auto& pt : ParticleTypeDatabase::DetectableTypes() ) {
        auto branch = ptype->AddBranch(pt);
        branch->AddHist1D(
                    [] (const ParticlePtr& p) { return p->Ek();},
                    HistogramFactory::Make1D("MC Energy " + pt->PrintName(),"E_{k} [MeV]","", energy_bins)
        );
        branch->AddHist2D(
                    [] (const ParticlePtr& p) { return make_tuple(p->Ek(), p->Theta()*TMath::RadToDeg());},
                    HistogramFactory::Make2D(
                "MC Energy " + pt->PrintName(),
                "E_{k} [MeV]",
                "#theta [#circ]",
                energy_bins,
                pt==&ParticleTypeDatabase::Proton ? theta_proton_bins : theta_bins)
        );
    }
}

ant::analysis::MCOverview::~MCOverview()
{

}

void ant::analysis::MCOverview::ProcessEvent(const ant::Event &event)
{
    const ParticleList& mc_particles = event.MCTrue().Particles().GetAll();
  //  cout << mc_particles.size() << endl;
    for( auto& mcp : mc_particles ) {
        mc_particle_stats.Fill(mcp);

    }
}

void ant::analysis::MCOverview::Finish()
{

}

void ant::analysis::MCOverview::ShowResult()
{
 /*
    canvas c("MC Overview");
    c << canvas::drawoption("colz");
    for(auto& plot : mc_particle_stats.Nodes()) {
        if( auto p = dynamic_pointer_cast<root_drawable_traits>(plot) ) {
            c << *p;
        }
    }
    c << canvas::cend;*/
}
