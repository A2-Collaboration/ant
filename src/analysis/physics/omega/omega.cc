#include "omega.h"
#include "data/Particle.h"
#include "data/Event.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include "plot/root_draw.h"
#include "plot/Histogram.h"
#include "utils/combinatorics.h"
#include "data/TaggerHit.h"
#include <string>
#include <iostream>
#include "plot/SmartHist.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;


void analysis::OmegaBase::ProcessEvent(const ant::Event &event)
{
    const auto& data = mode==DataMode::Reconstructed ? event.Reconstructed() : event.MCTrue();
    Analyse(data, event);
}

double OmegaBase::calcEnergySum(const ParticleList &particles) const
{
    double esum = 0.0;

    for( const ant::ParticlePtr& track : particles) {
        if( geo.DetectorFromAngles(track->Theta(), track->Phi()) == detector_t::NaI ) {
            esum += track->Ek();
        }
    }

    return esum;
}

ParticleList OmegaBase::getGeoAccepted(const ParticleList &p) const
{
    ParticleList list;
    for( auto& particle : p) {
        if( geo.DetectorFromAngles(particle->Theta(), particle->Phi()) != detector_t::None )
                list.emplace_back(particle);
    }
    return list;
}

OmegaBase::OmegaBase(const string &name, OmegaBase::DataMode m):
    Physics(name), mode(m)
{

}

void OmegaBase::Finish()
{
}

void OmegaBase::ShowResult()
{
}

//======= Omega Eta Gamma =====================================================================


OmegaEtaG::OmegaEtaG(OmegaBase::DataMode m):
    OmegaBase("OmegaEtaG_"+to_string(m), m)
{
    const ant::BinSettings imsacle(1000);

    ggg_gg     = HistFac.makeTH2D("3#gamma IM vs 2#gamma sub IM (signal only)","3#gamma IM [MeV]", "2#gamma sub IM [MeV]",imsacle,imsacle,"ggg_gg_omega");
    ggg_gg_bg  = HistFac.makeTH2D("3#gamma IM vs 2#gamma sub IM (background only)","3#gamma IM [MeV]", "2#gamma sub IM [MeV]",imsacle,imsacle,"ggg_gg_bg");
    ggg_gg_all  = HistFac.makeTH2D("3#gamma IM vs 2#gamma sub IM","3#gamma IM [MeV]", "2#gamma sub IM [MeV]",imsacle,imsacle,"ggg_gg_all");

    ggg = HistFac.makeTH1D("3#gamma IM","3#gamma IM [MeV]","",imsacle,"ggg");
    ggg_omega = HistFac.makeTH1D("3#gamma IM (from #omega)","3#gamma IM [MeV]","",imsacle,"ggg_omega");
    ggg_bg = HistFac.makeTH1D("3#gamma IM (non #omega)","3#gamma IM [MeV]","",imsacle,"ggg_bg");

    ggg_gg_omega_eta = HistFac.makeTH2D("3#gamma IM vs 2#gamma sub IM (#omega #rightarrow #eta #gamma only)","3#gamma IM [MeV]", "2#gamma sub IM [MeV]",imsacle,imsacle,"ggg_gg_omega_eta");
    ggg_gg_omega_pi0 = HistFac.makeTH2D("3#gamma IM vs 2#gamma sub IM (#omega #rightarrow #pi0 #gamma only)","3#gamma IM [MeV]", "2#gamma sub IM [MeV]",imsacle,imsacle,"ggg_gg_omega_pi0");

    ggg_omega_pi0oreta = HistFac.makeTH1D("3#gamma IM (#omega #rightarrow #pi^{0}/#eta)","3#gamma IM [MeV]","",imsacle,"ggg_omega_pi0oreta");

    steps = HistFac.makeTH1D("steps", "", "", BinSettings(10));
}

void OmegaEtaG::Analyse(const Event::Data &data, const Event &event)
{


    steps->Fill("Events seen",1);

    const auto nPhotons = data.Particles().Get(ParticleTypeDatabase::Photon).size();
    const auto nProtons = data.Particles().Get(ParticleTypeDatabase::Proton).size();

    if( nPhotons != 3 )    // require 3 photons
        return;
    steps->Fill("nPhotons",1);

    if( nProtons > 1)       // 0 or 1 proton
        return;

    steps->Fill("nProtons",1);

    if( data.Tracks().size() > 4 || data.Tracks().size() < 3 ) // not more then that (3photns +0or1 protons)
        return;

    steps->Fill("nTracks",1);

    const double CBESum = mode==DataMode::Reconstructed ? calcEnergySum(data.Particles().Get(ParticleTypeDatabase::Photon)) : data.TriggerInfos().CBEenergySum();

    if( CBESum < 550.0 )
        return;

    steps->Fill("ESum",1);

    const ParticleList photons = getGeoAccepted(data.Particles().Get(ParticleTypeDatabase::Photon));
    const ParticleList protons = getGeoAccepted(data.Particles().Get(ParticleTypeDatabase::Proton));

    auto& mctrue_final  = event.MCTrue().Particles().GetAll();
    auto& mctrue_interm = event.MCTrue().Intermediates().GetAll();

    bool is_omega_decay = false;
    const ParticleTypeDatabase::Type* subtype = nullptr;

    if(!mctrue_final.empty() && mctrue_final.at(0)->Type() == ParticleTypeDatabase::Proton) {

        if( mctrue_interm.size() >= 1) {
            is_omega_decay = (mctrue_interm.at(0)->Type() == ParticleTypeDatabase::Omega);

            if(mctrue_interm.size() ==2) {
            if( mctrue_interm.at(1)->Type() == ParticleTypeDatabase::Eta ||
                    mctrue_interm.at(1)->Type() == ParticleTypeDatabase::Pi0) {
                subtype = &(mctrue_interm.at(1)->Type());
            }
            }

        }
    }

    for( auto comb = makeCombination(photons,3); !comb.Done(); ++comb) {

        ParticleList ggg_list;
        ggg_list.assign(comb.begin(),comb.end());

        TLorentzVector gggState = *comb.at(0)+*comb.at(1)+*comb.at(2);
        const double gggIM = gggState.M();

        ggg->Fill(gggIM);

        if(is_omega_decay)
            ggg_omega->Fill(gggIM);
        else
            ggg_bg->Fill(gggIM);

        if(subtype!=nullptr) {
            ggg_omega_pi0oreta->Fill(gggIM);
        }

        for( auto gcomb = makeCombination(ggg_list,2); !gcomb.Done(); ++gcomb) {

            const TLorentzVector g1(*gcomb.at(0));
            const TLorentzVector g2(*gcomb.at(1));
            const TLorentzVector ggState = g1 + g2;
            const double ggIM = ggState.M();

            if(is_omega_decay) {
                ggg_gg->Fill(gggIM,ggIM);

                if(subtype == &ParticleTypeDatabase::Eta)
                    ggg_gg_omega_eta->Fill(gggIM,ggIM);
                else if(subtype == &ParticleTypeDatabase::Pi0)
                    ggg_gg_omega_pi0->Fill(gggIM,ggIM);
            }
            else
                ggg_gg_bg->Fill(gggIM,ggIM);

            ggg_gg_all->Fill(gggIM,ggIM);
        }

    }


}

void OmegaEtaG::ShowResult()
{
    canvas("OmegaEtaG Results")
            << ggg_gg << ggg_gg << ggg_gg_all
            << ggg
            << endc;
}


string to_string(const OmegaBase::DataMode &m)
{
    if(m == OmegaBase::DataMode::MCTrue) {
        return "MCTrue";
    } else {
        return "Reconstructed";
    }
}
