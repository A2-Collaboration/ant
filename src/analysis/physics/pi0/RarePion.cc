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
    const BinSettings im_binning(100,70,200);
    const BinSettings npart_binning(6,0,6);
    const BinSettings energy_binning(250,0,250);
    const BinSettings theta_binning(180,0,180);

    ggim     = HistFac.makeTH1D("2 #gamma IM",       "M_{#gamma #gamma} [MeV]",         "#", im_binning,    "IM_2g");
    gggim    = HistFac.makeTH1D("3 #gamma IM",       "M_{#gamma #gamma #gamma} [MeV]",  "#", im_binning,    "IM_3g");
    ggggim   = HistFac.makeTH1D("4 #gamma IM",       "M_{4#gamma} [MeV]",               "#", im_binning,    "IM_4g");
    im_3g    = HistFac.makeTH1D("IM of exact 3#gamma", "M_{3#gamma}",                   "#", im_binning,    "IM_3g_exact");

    nphotons = HistFac.makeTH1D("Number of photons", "N",                               "#", npart_binning, "N_g");
    nphotons_anglecut = HistFac.makeTH1D("Number of photons after cut on angle", "N",   "#", npart_binning, "N_g_cut");

    theta_vs_En = HistFac.makeTH2D("Theta vs Energy","E_{kin} [MeV]","#Theta [#circ]",energy_binning,theta_binning,"EnTheta");

    minAngle = HistFac.makeTH1D("angle between any 2 #gamma",       "#alpha [#circ]","#",theta_binning,"alpha");
    minAngle2g = HistFac.makeTH1D("angle between exact 2 #gamma",   "#alpha [#circ]","#",theta_binning,"alpha2g");
    minAngle3g = HistFac.makeTH1D("angle between 2 #gamma out of 3","#alpha [#circ]","#",theta_binning,"alpha3g");
    splitoffangle = HistFac.makeTH1D("smallest angle between 2 #gamma","#alpha [#circ]","#",theta_binning,"alphamin");


    // Build a map of ParticleType -> Histogram, and fill it
//    for( auto& type : ParticleTypeDatabase() ) {
//        EHists[&type] = HistFac.makeTH1D(type.PrintName()+" Energy", "E [MeV]", "", energy_binning);
//    }

}


void RarePion::ProcessEvent(const TEvent& event, manager_t&)

{
    const TParticleList& photons = event.Reconstructed->Particles.Get(ParticleTypeDatabase::Photon);
    const TParticleList& all = event.Reconstructed->Particles.GetAll();

//    for( auto& particle : all ) {

//        // fill the histogram corresponding to the partice type of the current particle
//        auto entry = EHists.find(&particle->Type());
//        if( entry != EHists.end()) {
//            entry->second->Fill(particle->Ek());
//        }

//    }

    nphotons->Fill(photons.size());

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

    for( auto& p: photons ) {
        theta_vs_En->Fill(p->Ek(),p->Theta()*TMath::RadToDeg());
    }

    for( auto comb = utils::makeCombination(photons,2); !comb.Done(); ++comb) {
        minAngle->Fill(comb.at(0)->Angle(comb.at(1)->Vect())*TMath::RadToDeg());
    }

    if(photons.size() == 2) { // cut already on # of clusters
        for( auto comb = utils::makeCombination(photons,2); !comb.Done(); ++comb) {
            minAngle2g->Fill(comb.at(0)->Angle(comb.at(1)->Vect())*TMath::RadToDeg());
        }
    }

    if(photons.size() == 3) { // cut already on # of clusters
        const TParticle pi0 ( ParticleTypeDatabase::Pi0, *photons.at(0) + *photons.at(1) + *photons.at(2));
        im_3g->Fill(pi0.M());
        Double_t anglemin = 180;
        for( auto comb = utils::makeCombination(photons,2); !comb.Done(); ++comb) {
            Double_t angle = comb.at(0)->Angle(comb.at(1)->Vect())*TMath::RadToDeg();
            //cout << "angle: " << angle << "\n";
            minAngle3g->Fill(angle);
            if(anglemin > angle) {
                anglemin = angle;
            }
        }
        splitoffangle->Fill(anglemin);
        if (anglemin > 30) {
            nphotons_anglecut->Fill(photons.size());
        }
    }

}

void RarePion::Finish()
{

}

void RarePion::ShowResult()
{
//    canvas cc("RarePion");
//    cc <<
//          ggim <<
//          gggim <<
//          ggggim <<
//          nphotons <<
//          //nprotons <<
//          endc;

//    canvas cc2("RarePion2");
//    for( auto& e : EHists) {
//        cc2 << e.second;
//    }
//    cc2 << endc;


}

AUTO_REGISTER_PHYSICS(RarePion)
