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
#include <algorithm>
#include <iostream>

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace ant::analysis::data;

void OmegaBase::ProcessEvent(const Event &event)
{
    const auto& data = mode==DataMode::Reconstructed ? event.Reconstructed() : event.MCTrue();
    Analyse(data, event);
}

double OmegaBase::calcEnergySum(const ParticleList &particles) const
{
    double esum = 0.0;

    for( const ParticlePtr& p : particles) {
        if( geo.DetectorFromAngles(p->Theta(), p->Phi()) == Detector_t::Type_t::CB ) {
            esum += p->Ek();
        }
    }

    return esum;
}

ParticleList OmegaBase::getGeoAccepted(const ParticleList &p) const
{
    ParticleList list;
    for( auto& particle : p) {
        if( geo.DetectorFromAngles(particle->Theta(), particle->Phi()) != Detector_t::Any_t::None )
                list.emplace_back(particle);
    }
    return list;
}

string OmegaBase::GetDecayString(const ParticleList &particles)
{
    stringstream s;
    if(! particles.empty()) {
        const auto& start = particles.front();
        Particle::RecPrint(start, s);
    }

    return s.str();
}

OmegaBase::OmegaBase(const string &name, PhysOptPtr opts):
    Physics(name, opts), mode(DataMode::Reconstructed)
{

}

void OmegaBase::Finish()
{
}

void OmegaBase::ShowResult()
{
}

//======= Omega Eta Gamma =====================================================================


OmegaEtaG::OmegaEtaG(PhysOptPtr opts):
    OmegaBase("OmegaEtaG", opts)
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

    steps = HistFac.makeTH1D("steps", "", "", BinSettings(10));

    if(GetOption("OmegaEtaGMode") == "McTrue") {
        mode = OmegaBase::DataMode::MCTrue;
    }

    VLOG(8) << "mode option is " << GetOption("OmegaEtaGMode");
}

OmegaEtaG::perDecayhists_t OmegaEtaG::makePerDecayHists(const string &title)
{
    perDecayhists_t h;

    auto pref(title);
    std::replace( pref.begin(), pref.end(), '[', '_');
    std::replace( pref.begin(), pref.end(), ']', '_');
    std::replace( pref.begin(), pref.end(), ' ', '_');
    std::replace( pref.begin(), pref.end(), '#', '_');

    h.gg = HistFac.makeTH1D("2#gamma "+title,"2#gamma IM [MeV]","#",imbinning,pref+"gg");
    h.ggg  = HistFac.makeTH1D("3#gamma "+title,"3#gamma IM [MeV]","",imbinning,pref+"ggg");
    h.mm = HistFac.makeTH1D("MM "+title,"MM [MeV]","",mmbinning,pref+"mm");
    h.angle_p = HistFac.makeTH1D(title+"Angle p rec/true","angle [#circ]","",BinSettings(4*360,0,180),pref+"angle_p");
    h.angle_p_ggg = HistFac.makeTH1D(title+"Angle p ggg","angle [#circ]","",BinSettings(4*360,0,180),pref+"angle_p_ggg");
    h.p_phi_diff = HistFac.makeTH1D(title+"p phi diff","angle [#circ]","",BinSettings(4*360,-90,90),pref+"p_phi_diff");

    return h;
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

    if( data.Candidates().size() > 4 || data.Candidates().size() < 3 ) // not more then that (3photns +0or1 protons)
        return;

    steps->Fill("nCandidates",1);

    const double CBESum = mode==DataMode::Reconstructed ? calcEnergySum(data.Particles().Get(ParticleTypeDatabase::Photon)) : data.TriggerInfos().CBEenergySum();

    if( CBESum < 550.0 )
        return;

    steps->Fill("ESum",1);

    const ParticleList photons = getGeoAccepted(data.Particles().Get(ParticleTypeDatabase::Photon));
    const ParticleList protons = getGeoAccepted(data.Particles().Get(ParticleTypeDatabase::Proton));


    bool is_omega_decay = false;
    const ParticleTypeDatabase::Type* subtype = nullptr;

    perDecayhists_t* h = nullptr;

    const string decaystring = GetDecayString(event.MCTrue().Intermediates().GetAll());

    for( auto comb = utils::makeCombination(photons,3); !comb.Done(); ++comb) {

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
                VLOG(9) << "Adding Histograms for " << decaystring;
                gg_decays[decaystring] = makePerDecayHists(decaystring);
                h = &(gg_decays[decaystring]);
            } else {
                h = &(entry->second);
            }
        }

        if(h)
            h->ggg->Fill(gggIM);

        const auto& mc_protons = event.MCTrue().Particles().Get(ParticleTypeDatabase::Proton);

        const ParticlePtr mc_p = !mc_protons.empty() ? mc_protons.at(0) : shared_ptr<Particle>(nullptr);

        if(h) {
            for(auto& th : event.MCTrue().TaggerHits()) {
                const TLorentzVector beam_target = th->PhotonBeam() + TLorentzVector(0, 0, 0, ParticleTypeDatabase::Proton.Mass());
                TLorentzVector mm = beam_target - gggState;
                mm.Boost(-beam_target.BoostVector());
                h->mm->Fill(mm.M());

                if(mc_p) {
                    TLorentzVector mc_p_v = TLorentzVector(*mc_p);
                    mc_p_v.Boost(-beam_target.BoostVector());
                    const auto angle = mc_p_v.Angle(mm.Vect());
                    h->angle_p->Fill(angle * TMath::RadToDeg());

                    TLorentzVector ggg_boost = gggState;
                    ggg_boost.Boost(-beam_target.BoostVector());
                    const auto angle_pggg = ggg_boost.Angle(mc_p_v.Vect());
                    h->angle_p_ggg->Fill(angle_pggg * TMath::RadToDeg());

                    const auto p_phi_diff = TVector2::Phi_mpi_pi(mm.Phi() - mc_p_v.Phi());
                    h->p_phi_diff->Fill(p_phi_diff);

                }
            }
        }

        for( auto gcomb = utils::makeCombination(ggg_list,2); !gcomb.Done(); ++gcomb) {

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

            if(h)
                h->gg->Fill(ggIM);

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

    std::list<perDecayhists_t*> histlist;

    for(auto& hist : gg_decays) {

       histlist.emplace_back(&(hist.second));
    }

    histlist.sort( [] (const perDecayhists_t* h1, const perDecayhists_t* h2) -> bool {return h1->gg->GetEntries() > h2->gg->GetEntries();});

    std::vector<Color_t> color = {kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kPink+9, kSpring+10, kGray};
    auto cit = circit<decltype(color.begin())>(color.begin(),color.end());

    hstack stack("decays","gg");
    hstack stack2("ggg","ggg");
    hstack stack3("mm","mm");
    hstack stack4("angle_p","angle_p");
    hstack stack5("angle_pggg","angle_pggg");
    hstack stack6("p_phi_diff","p_phi_diff");

    for(auto& hist:histlist) {
        hist->gg->SetFillColor(*cit);
        hist->ggg->SetFillColor(*cit);
        hist->mm->SetFillColor(*cit);
        cit.next();
        stack << hist->gg;
        stack2 << hist->ggg;
        stack3 << hist->mm;
        stack4 << hist->angle_p;
        stack5 << hist->angle_p_ggg;
        stack6 << hist->p_phi_diff;
    }

    canvas("OmegaEtaG per Decay Results gg") << stack << endc;
    canvas("OmegaEtaG per Decay Results ggg") << stack2 << endc;
    canvas("OmegaEtaG per Decay Results mm") << drawoption("pads") << stack3 << endc;
    canvas("OmegaEtaG per Decay Results Angle p") << drawoption("pads") << stack4 << endc;
    canvas("OmegaEtaG per Decay Results Angle p ggg") << drawoption("pads") << stack5 << endc;
    canvas("OmegaEtaG per Decay Resultsp phi diff") << drawoption("pads") << stack6 << endc;
}


string to_string(const OmegaBase::DataMode &m)
{
    if(m == OmegaBase::DataMode::MCTrue) {
        return "MCTrue";
    } else {
        return "Reconstructed";
    }
}


AUTO_REGISTER_PHYSICS(OmegaEtaG, "OmegaEtaG")
