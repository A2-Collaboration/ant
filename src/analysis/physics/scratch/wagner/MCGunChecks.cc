#include "MCGunChecks.h"

#include "utils/ParticleTools.h"
#include "tree/TParticle.h"

#include "plot/HistStyle.h"

#include "base/Logger.h"
#include "base/std_ext/math.h"

#include "TH1D.h"
#include "TTree.h"

#include "utils/Combinatorics.h"

#include <cmath>
#include <iostream>

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;


scratch_wagner_MCGunChecks::scratch_wagner_MCGunChecks(const std::string& name, OptionsPtr opts):
    Physics(name,opts)
{
    t.CreateBranches(HistFac.makeTTree("tree"));
}


void scratch_wagner_MCGunChecks::ProcessEvent(const TEvent& event, manager_t&)
{
    if (!event.MCTrue().ParticleTree)
        return;

    const auto particles = [](const TEvent& event)
    {
        const auto gPT = event.MCTrue().ParticleTree->Daughters();
        vector<TParticlePtr> ret(gPT.size());
        transform(gPT.begin(), gPT.end(), ret.begin(),
                  [](const TParticleTree_t& en){ return en->Get(); }
        );
        return ret;
    }(event);

    for (const auto& gP: particles)
    {
        t.names().emplace_back(gP->Type().Name());
        t.thetas_true().push_back(gP->Theta());
        t.phis_true().push_back(gP->Phi());
    }
    t.multiplicity = event.Reconstructed().Candidates.size();

    if (t.multiplicity == 1) {
        auto& cand = event.Reconstructed().Candidates.front();
        t.energies().push_back(cand.CaloEnergy);
        t.thetas().push_back(cand.Theta);
        t.phis().push_back(cand.Phi);
        auto particle_list = utils::ParticleTypeList::Make(event.MCTrue().ParticleTree).GetAll();

        if (particle_list.size() != 1)
            throw runtime_error("Too many true particles");
        t.energies_true().push_back(particle_list.front()->Ek());
        //t.thetas_true().push_back(particle_list.front()->Theta());
        //t.phis_true().push_back(particle_list.front()->Phi());
    }


    t.fillAndReset();

}


void scratch_wagner_MCGunChecks::Finish()
{

}

void scratch_wagner_MCGunChecks::ShowResult()
{
    auto tree = t.Tree;
    const auto e_sigma = BinSettings(600, -150, 150);
    canvas(string("Energy Resolution ") + t.names().front())
            << TTree_drawable(tree, "multiplicity")
            << TTree_drawable(tree, "energies - energies_true", "",  // no cuts
                              "Energy Resolution MC Gun", "#sigma(E)", "#", e_sigma)
            << drawoption("colz")
            << TTree_drawable(tree, "energies - energies_true:thetas*180/3.1415",
                              "(energies - energies_true) > -150 && (energies - energies_true) < 150")
            << TTree_drawable(tree, "energies - energies_true:energies_true",
                              "(energies - energies_true) > -150 && (energies - energies_true) < 150")
            << endc;

    canvas(string("Theta Resolution ") + t.names().front())
            << TTree_drawable(tree, "multiplicity")
            << TTree_drawable(tree, "(thetas - thetas_true)*180/3.1415", "",
                              "Theta Resolution MC Gun", "#sigma(#vartheta)", "#", BinSettings(100, -25, 25))
            << drawoption("colz")
            << TTree_drawable(tree, "(thetas - thetas_true)*180/3.1415:energies",
                              "(thetas - thetas_true)*180/3.1415 > -25 && (thetas - thetas_true)*180/3.1415 < 25")
            << TTree_drawable(tree, "(thetas - thetas_true)*180/3.1415:thetas_true*180/3.1415",
                              "(thetas - thetas_true)*180/3.1415 > -25 && (thetas - thetas_true)*180/3.1415 < 25")
            << endc;
}


AUTO_REGISTER_PHYSICS(scratch_wagner_MCGunChecks)
