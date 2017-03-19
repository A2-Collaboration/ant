#include "RarePion.h"

#include "base/ParticleType.h"
#include "utils/Combinatorics.h"
#include "utils/ParticleTools.h"

#include "TCanvas.h"
#include "TH1D.h"

#include <iostream>
#include <memory>

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
    im_2g    = HistFac.makeTH1D("IM of exact 2#gamma", "M_{2#gamma}",                   "#", im_binning,    "IM_2g_exact");

    nphotons = HistFac.makeTH1D("Number of photons", "N",                               "#", npart_binning, "N_g");
    nphotons_anglecut = HistFac.makeTH1D("Number of photons after cut on angle > 40#circ", "N",   "#", npart_binning, "N_g_cut");

    theta_vs_En = HistFac.makeTH2D("Theta vs Energy","E_{kin} [MeV]","#Theta [#circ]",energy_binning,theta_binning,"EnTheta");
    splitoffenergy = HistFac.makeTH2D("Energy of 2 nearby photons (out of exact 3)","E_{kin} of one #gamma [MeV]","E_{kin} of the other #gamma [MeV]",energy_binning,energy_binning,"EnSplit");

    minAngle = HistFac.makeTH1D("angle between any 2 #gamma of 3 clusters",       "#alpha [#circ]","#",theta_binning,"alpha");
    minAngle2g = HistFac.makeTH1D("angle between exact 2 #gamma (2 clusters)",   "#alpha [#circ]","#",theta_binning,"alpha2g");
    splitoffangle = HistFac.makeTH1D("smallest angle between 2 #gamma (3 clusters)","#alpha [#circ]","#",theta_binning,"alphamin");

    // Build a map of ParticleType -> Histogram, and fill it
//    for( auto& type : ParticleTypeDatabase() ) {
//        EHists[&type] = HistFac.makeTH1D(type.PrintName()+" Energy", "E [MeV]", "", energy_binning);
//    }

}


void RarePion::ProcessEvent(const TEvent& event, manager_t&)
{
    auto recon_particles = utils::ParticleTypeList::Make(event.Reconstructed().Candidates);
    const TParticleList& photons = recon_particles.Get(ParticleTypeDatabase::Photon);
    //const TParticleList& all = event.Reconstructed().Particles.GetAll();

//    for( auto& particle : all ) {

//        // fill the histogram corresponding to the partice type of the current particle
//        auto entry = EHists.find(&particle->Type());
//        if( entry != EHists.end()) {
//            entry->second->Fill(particle->Ek());
//        }

//    }

    constexpr auto radtodeg = std_ext::radian_to_degree(1.0);

    nphotons->Fill(photons.size());

    Bool_t checkCB = true;
    for( auto& p: photons ) {
        if (! (p->Candidate->Detector & Detector_t::Type_t::CB)){
            checkCB = false;
        }
    }

    if ( checkCB == true ){ // only continue when all photons are in CB
        for( auto& p: photons ) {
            theta_vs_En->Fill(p->Ek(),p->Theta()*radtodeg);
        }

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

        auto combinations4 = utils::makeCombination(photons,4);
        do {
            LorentzVec v;
            for( auto& i: combinations4 ) {
                v += *i;
            }

            ggggim->Fill(v.M());

        } while(combinations4.next());

        for( auto comb = utils::makeCombination(photons,2); !comb.done(); ++comb) {
            minAngle->Fill(comb.at(0)->Angle(comb.at(1)->p)*radtodeg);
        }

        if(photons.size() == 2) { // cut already on # of clusters
            const TParticle pi02g ( ParticleTypeDatabase::Pi0, *photons.at(0) + *photons.at(1));
            im_2g->Fill(pi02g.M());
            for( auto comb = utils::makeCombination(photons,2); !comb.done(); ++comb) {
                minAngle2g->Fill(comb.at(0)->Angle(comb.at(1)->p)*radtodeg);
            }
        }

        if(photons.size() == 3) { // cut already on # of clusters
            const TParticle pi03g ( ParticleTypeDatabase::Pi0, *photons.at(0) + *photons.at(1) + *photons.at(2));
            im_3g->Fill(pi03g.M());
            Double_t anglemin = 180;
            angle01 = photons.at(0)->Angle(photons.at(1)->p)*radtodeg;
            angle02 = photons.at(0)->Angle(photons.at(2)->p)*radtodeg;
            angle12 = photons.at(1)->Angle(photons.at(2)->p)*radtodeg;
            if (angle01 < angle02){
                if (angle01 < angle12){
                    //Fall1, vectors 01
                    splitoffenergy->Fill(photons.at(0)->Ek(),photons.at(1)->Ek());
                    anglemin = angle01;
                }
                else {
                    //Fall2, vectors 12
                    splitoffenergy->Fill(photons.at(1)->Ek(),photons.at(2)->Ek());
                    anglemin = angle12;
                }
            }
            else  if (angle02 < angle12){
                //Fall3, vectors 02
                splitoffenergy->Fill(photons.at(0)->Ek(),photons.at(2)->Ek());
                anglemin = angle02;
            }
            splitoffangle->Fill(anglemin);
            if (anglemin > 40) {
                nphotons_anglecut->Fill(photons.size());
            }
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
