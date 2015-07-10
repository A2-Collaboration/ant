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
#include <sstream>
#include "TH3.h"
#include "base/Logger.h"

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

string OmegaBase::GetDecayString(const ParticleList &particles)
{
    stringstream s;
    for(auto& p : particles) {
        s << p->Type().PrintName() << " ";
    }
    return s.str();
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
    ggg_gg     = HistFac.makeTH2D("3#gamma IM vs 2#gamma sub IM (signal only)","3#gamma IM [MeV]", "2#gamma sub IM [MeV]",imbinning,imbinning,"ggg_gg_omega");
    ggg_gg_bg  = HistFac.makeTH2D("3#gamma IM vs 2#gamma sub IM (background only)","3#gamma IM [MeV]", "2#gamma sub IM [MeV]",imbinning,imbinning,"ggg_gg_bg");
    ggg_gg_all  = HistFac.makeTH2D("3#gamma IM vs 2#gamma sub IM","3#gamma IM [MeV]", "2#gamma sub IM [MeV]",imbinning,imbinning,"ggg_gg_all");

    ggg = HistFac.makeTH1D("3#gamma IM","3#gamma IM [MeV]","",imbinning,"ggg");
    ggg_omega = HistFac.makeTH1D("3#gamma IM (from #omega)","3#gamma IM [MeV]","",imbinning,"ggg_omega");
    ggg_bg = HistFac.makeTH1D("3#gamma IM (non #omega)","3#gamma IM [MeV]","",imbinning,"ggg_bg");

    ggg_gg_omega_eta = HistFac.makeTH2D("3#gamma IM vs 2#gamma sub IM (#omega #rightarrow #eta #gamma only)","3#gamma IM [MeV]", "2#gamma sub IM [MeV]",imbinning,imbinning,"ggg_gg_omega_eta");
    ggg_gg_omega_pi0 = HistFac.makeTH2D("3#gamma IM vs 2#gamma sub IM (#omega #rightarrow #pi0 #gamma only)","3#gamma IM [MeV]", "2#gamma sub IM [MeV]",imbinning,imbinning,"ggg_gg_omega_pi0");

    ggg_omega_pi0oreta = HistFac.makeTH1D("3#gamma IM (#omega #rightarrow #pi^{0}/#eta)","3#gamma IM [MeV]","",imbinning,"ggg_omega_pi0oreta");

    ggg_gg_decays = HistFac.makeTH3D(
                "3#gamma IM vs 2#gamma sub IM (signal only)","3#gamma IM [MeV]", "2#gamma sub IM [MeV]","decay",
                imbinning,imbinning,BinSettings(1),"ggg_gg_decays");

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

    TH1D* h = nullptr;

    const string decaystring = GetDecayString(event.MCTrue().Intermediates().GetAll());

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

        if(omega_range.Contains(gggIM) && h == nullptr) {

            auto entry = gg_decays.find(decaystring);
            if(entry == gg_decays.end()) {
                VLOG(9) << "Adding Histogram for " << decaystring;
                h = HistFac.makeTH1D(decaystring,"2#gamma IM [MeV]","#",imbinning);
                gg_decays[decaystring] = h;
            } else {
                h = entry->second;
            }
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

            ggg_gg_decays->Fill(gggIM,ggIM,decaystring.c_str(),1);

            if(h)
                h->Fill(ggIM);

        }

    }


}

template <typename T>
class circit{
protected:
    T s;
    T e;
    T p;
public:
    circit(const T& start, const T& stop):s(start),e(stop),p(start) {}

    void next() {
        ++p;
        if(p==e)
            p=s;
    }

    void previous() {
        if(p==s) {
            p=e;
        }
        --p;
    }

    void Reset() {
        p = s;
    }

    typename T::value_type operator *() { return *p; }

};

void OmegaEtaG::ShowResult()
{
    canvas("OmegaEtaG Results")
            << ggg_gg << ggg_gg << ggg_gg_all
            << ggg
            << endc;

    canvas c1("OmegaEtaG per Decay Results");
    hstack stack("decays","Decays");

    std::list<TH1D*> histlist;

    for(auto& hist : gg_decays) {

       histlist.emplace_back(hist.second);
    }

    histlist.sort( [] (const TH1D* h1, const TH1D* h2) -> bool {return h1->GetEntries() > h2->GetEntries();});
    std::vector<Color_t> color = {kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kPink+9, kSpring+10, kGray};

    auto cit = circit<decltype(color.begin())>(color.begin(),color.end());

    for(auto& hist:histlist) {
        hist->SetFillColor(*cit);
        cit.next();
        stack << hist;
    }

    c1 << stack << endc;

}


string to_string(const OmegaBase::DataMode &m)
{
    if(m == OmegaBase::DataMode::MCTrue) {
        return "MCTrue";
    } else {
        return "Reconstructed";
    }
}
