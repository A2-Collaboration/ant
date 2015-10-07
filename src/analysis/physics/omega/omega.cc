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
#include "TH3.h"
#include "base/Logger.h"
#include <algorithm>
#include <iostream>
#include "base/std_ext/math.h"

#include "TTree.h"
#include "base/iterators.h"

#include "utils/particle_tools.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace ant::analysis::data;
using namespace ant::std_ext;

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
    h.calc_proton_energy_theta = HistFac.makeTH2D(title+" Calc Proton","E [MeV]","#Theta [#circ]",BinSettings(1000),BinSettings(360,0,180),pref+"calc_p");
    h.calc_proton_special = HistFac.makeTH2D(title+" Calc Proton Special","E [MeV]","#Theta [#circ]",BinSettings(1000),BinSettings(360,0,180),pref+"calc_p_special");
    h.nCand = HistFac.makeTH1D(title+": # Photons","Photons/Event","",BinSettings(15),pref+"ncand");
    return h;
}

TLorentzVector Boost(const TLorentzVector& lv, const TVector3& boost) {
    TLorentzVector b(lv);
    b.Boost(boost);
    return b;
}

void OmegaEtaG::Analyse(const Event::Data &data, const Event &event)
{

    steps->Fill("Events seen",1);

    const auto nPhotons = data.Particles().Get(ParticleTypeDatabase::Photon).size();
    const auto nProtons = data.Particles().Get(ParticleTypeDatabase::Proton).size();

    perDecayhists_t* h = nullptr;

    const string decaystring = utils::ParticleTools::GetDecayString(event.MCTrue().Intermediates().GetAll());

    if( h == nullptr) {

        auto entry = gg_decays.find(decaystring);
        if(entry == gg_decays.end()) {
            VLOG(9) << "Adding Histograms for " << decaystring;
            gg_decays[decaystring] = makePerDecayHists(decaystring);
            h = &(gg_decays[decaystring]);
        } else {
            h = &(entry->second);
        }
    }

    if(h) {
        h->nCand->Fill(nPhotons);
    }

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



        if(h) {
            h->ggg->Fill(gggIM);
        }

        const auto& mc_protons = event.MCTrue().Particles().Get(ParticleTypeDatabase::Proton);

        const ParticlePtr mc_p = !mc_protons.empty() ? mc_protons.at(0) : shared_ptr<Particle>(nullptr);

        if(omega_range.Contains(gggIM) && h) {
            for(auto& th : event.MCTrue().TaggerHits()) {
                const TLorentzVector beam_target = th->PhotonBeam() + TLorentzVector(0, 0, 0, ParticleTypeDatabase::Proton.Mass());
                const TLorentzVector mm = beam_target - gggState;
                const TLorentzVector mm_boosted = Boost(mm,-beam_target.BoostVector());


                h->calc_proton_energy_theta->Fill(mm.Energy()-ParticleTypeDatabase::Proton.Mass(), mm.Theta() * TMath::RadToDeg());

                h->mm->Fill(mm.M());

                if(mc_p) {

                    const Particle p_special(ParticleTypeDatabase::Proton, mm.Energy()-ParticleTypeDatabase::Proton.Mass(), mc_p->Theta(), mc_p->Phi());
                    h->calc_proton_special->Fill(p_special.Ek(), p_special.Theta()*TMath::RadToDeg());

                    TLorentzVector mc_p_v = TLorentzVector(*mc_p);
                    mc_p_v.Boost(-beam_target.BoostVector());
                    const auto angle = mc_p_v.Angle(mm_boosted.Vect());
                    h->angle_p->Fill(angle * TMath::RadToDeg());

                    TLorentzVector ggg_boost = gggState;
                    ggg_boost.Boost(-beam_target.BoostVector());
                    const auto angle_pggg = ggg_boost.Angle(mc_p_v.Vect());
                    h->angle_p_ggg->Fill(angle_pggg * TMath::RadToDeg());

                    const auto p_phi_diff = TVector2::Phi_mpi_pi(mm_boosted.Phi() - mc_p_v.Phi());
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

            if(h && omega_range.Contains(gggIM))
                h->gg->Fill(ggIM);

        }

    }


}

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


    auto cit = getCirculatIterator(ColorPalette::Colors.begin(), ColorPalette::Colors.end());

    hstack stack("decays","gg");
    hstack stack2("ggg","ggg");
    hstack stack3("mm","mm");
    hstack stack4("angle_p","angle_p");
    hstack stack5("angle_pggg","angle_pggg");
    hstack stack6("p_phi_diff","p_phi_diff");
    canvas c7("OmegaEtaG per decay Results calc p");
    c7 << drawoption("colz");
    canvas c8("OmegaEtaG per decay Results calc p special");
    c8 << drawoption("colz");

    int i=0;
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
        c7 << hist->calc_proton_energy_theta;
        c8 << hist->calc_proton_special;
        if(i++==9)
            break;
    }

    canvas("OmegaEtaG per Decay Results gg") << stack << endc;
    canvas("OmegaEtaG per Decay Results ggg") << stack2 << endc;
    canvas("OmegaEtaG per Decay Results mm") << drawoption("pads") << stack3 << endc;
    canvas("OmegaEtaG per Decay Results Angle p") << drawoption("pads") << stack4 << endc;
    canvas("OmegaEtaG per Decay Results Angle p ggg") << drawoption("pads") << stack5 << endc;
    canvas("OmegaEtaG per Decay Resultsp phi diff") << drawoption("pads") << stack6 << endc;
    c7 << endc;
    c8 << endc;
}


string to_string(const OmegaBase::DataMode &m)
{
    if(m == OmegaBase::DataMode::MCTrue) {
        return "MCTrue";
    } else {
        return "Reconstructed";
    }
}




OmegaMCTruePlots::PerChannel_t::PerChannel_t(const string& Title, SmartHistFactory& hf):
    title(Title)
{
    proton_E_theta = hf.makeTH2D(title,"E [MeV]","#theta [#circ]",BinSettings(1000),BinSettings(360,0,180), title+"_e_theta");
}

void OmegaMCTruePlots::PerChannel_t::Show()
{
    canvas("Omega per Channel: "+title) << drawoption("colz") << proton_E_theta << endc;
}

void OmegaMCTruePlots::PerChannel_t::Fill(const Event::Data& d)
{
    const auto& protons = d.Particles().Get(ParticleTypeDatabase::Proton);
    if(!protons.empty()) {
        const auto& p = protons.at(0);
        proton_E_theta->Fill(p->Ek(), p->Theta()*TMath::RadToDeg());
    }
}



OmegaMCTruePlots::OmegaMCTruePlots(PhysOptPtr opts):
    Physics("OmegaMCTruePlots", opts)
{

}

void OmegaMCTruePlots::ProcessEvent(const Event& event)
{
    const auto& decaystring =utils::ParticleTools::GetProductionChannelString(event.MCTrue().Intermediates().GetAll());

    auto e = channels.find(decaystring);

    if(e == channels.end()) {
        channels.insert({decaystring, PerChannel_t(decaystring,HistFac)});
    }

    e = channels.find(decaystring);
    e->second.Fill(event.MCTrue());
}

void OmegaMCTruePlots::Finish()
{

}

void OmegaMCTruePlots::ShowResult()
{
    canvas c("OmegaMCTrue p E Theta");
    c << drawoption("colz");

    list<TH2D*> hists;
    for(auto& entry : channels) {
        hists.push_back(entry.second.proton_E_theta);
    }

    hists.sort([](const TH2D* a, const TH2D* b) {return a->GetEntries() > b->GetEntries();});

    int i=0;
    for(auto& h : hists) {
        c << h;
        i++;
        if(i>=9)
            break;
    }

    c << endc;
}



TLorentzVector OmegeMCTree::getGamma1() const
{
    return gamma1_vector;
}

void OmegeMCTree::setGamma1(const TLorentzVector& value)
{
    gamma1_vector = value;
}

OmegeMCTree::OmegeMCTree(PhysOptPtr opts): Physics("OmegaMCTree", opts) {
    tree=new TTree("omegatree","omgega eta gamma MC true");
    tree->Branch("p", &proton_vector);
    tree->Branch("omega", &omega_vector);
    tree->Branch("gamma1", &gamma1_vector);
    tree->Branch("eta", &eta_vector);
    tree->Branch("gamma2", &gamma2_vector);
    tree->Branch("gamma3", &gamma3_vector);
}

OmegeMCTree::~OmegeMCTree()
{

}

void OmegeMCTree::ProcessEvent(const Event& event)
{
    proton_vector.SetPxPyPzE(0,0,0,0);
    omega_vector.SetPxPyPzE(0,0,0,0);
    gamma1_vector.SetPxPyPzE(0,0,0,0);
    eta_vector.SetPxPyPzE(0,0,0,0);
    gamma2_vector.SetPxPyPzE(0,0,0,0);
    gamma3_vector.SetPxPyPzE(0,0,0,0);

    const auto& bpl = event.MCTrue().Intermediates().Get(ParticleTypeDatabase::BeamProton);
    if(bpl.size() == 1) {
        const ParticlePtr& bp = bpl.at(0);

        if(bp->Daughters().size() == 2) {
            for(const auto& d : bp->Daughters()) {
                if(d->Type() == ParticleTypeDatabase::Proton) {
                    proton_vector = *d;
                } else if(d->Type() == ParticleTypeDatabase::Omega) {
                    omega_vector = *d;
                    if(d->Daughters().size() ==2 ) {
                        for(const ParticlePtr& e : d->Daughters()) {
                            if(e->Type() == ParticleTypeDatabase::Eta) {
                                eta_vector = *e;
                                if(e->Daughters().size() == 2) {
                                    for(const ParticlePtr& f : e->Daughters()) {
                                        if(f->Type() == ParticleTypeDatabase::Photon) {
                                            if(gamma2_vector.E() ==0) {
                                                gamma2_vector = *f;
                                            } else {
                                                if(gamma3_vector.E()==0) {
                                                    gamma3_vector = *f;
                                                }
                                            }
                                        }
                                    }
                                }
                            } else if(e->Type() == ParticleTypeDatabase::Photon) {
                                gamma1_vector = *e;
                            }
                        }
                    }
                }
            }
        }

        if(omega_vector.E() != 0 && proton_vector.E() !=0 && gamma1_vector.E() !=0 && eta_vector.E() !=0 && gamma2_vector.E() !=0 && gamma3_vector.E() != 0) {
            tree->Fill();
        } else {
            LOG(WARNING) << "not complete";
        }
    }
}

void OmegeMCTree::ShowResult()
{

}


template <typename it_type>
TLorentzVector LVSum(it_type begin, it_type end) {
    TLorentzVector v;

    while(begin!=end) {
        v += **begin;
        ++begin;
    }

    return v;
}

template <typename it_type>
double TimeAvg(it_type begin, it_type end) {
    double v=0.0;

    while(begin!=end) {
        const ParticlePtr& p = *begin;
        v += p->Candidates().at(0)->Time();
        ++begin;
    }

    return v;
}

void OmegaEtaG2::Analyse(const Event::Data &data, const Event &event)
{
    const ParticleList& iphotons = data.Particles().Get(ParticleTypeDatabase::Photon);
    const ParticleList& iprotons = (data_proton ? event.Reconstructed() : event.MCTrue()).Particles().Get(ParticleTypeDatabase::Proton);

    if(iphotons.size() != 3)
        return;

    if(iprotons.size() != 1)
        return;

    const ParticleList protons = getGeoAccepted(iprotons);

    if(protons.size() != 1)
        return;

    const ParticleList photons = getGeoAccepted(iphotons);

    if(photons.size() != 3)
        return;

    const auto Esum = calcEnergySum2(data);

    if(Esum < 550.0)
        return;

    const auto& proton = protons.at(0);

    pEk    = proton->Ek();
    pTheta = radian_to_degree(proton->Theta());
    pPhi   = radian_to_degree(proton->Phi());
    pTime  = data_proton ? proton->Candidates().at(0)->Time() : 0.0;

    const TLorentzVector ggg = LVSum(photons.begin(), photons.end());
    gggTime  = TimeAvg(photons.begin(), photons.end());
    gggIM    = ggg.M();
    gggTheta = radian_to_degree(ggg.Theta());
    gggPhi   = radian_to_degree(ggg.Phi());

    TLorentzVector gg;
    int ngg=0;
    for( auto comb = utils::makeCombination(photons,2); !comb.Done(); ++comb) {
        gg.SetPtEtaPhiE(0,0,0,0);
        auto i = comb.begin();
        gg+=**i;
        ++i;
        gg+=**i;

        ggIM[ngg++] = gg.M();
    }

    rf = -1; // TODO

    for(const TaggerHitPtr& t : data_tagger?data.TaggerHits():event.MCTrue().TaggerHits()) {
        tagch   = t->Channel();
        tagtime = t->Time();

        const TLorentzVector beam_target = t->PhotonBeam() + TLorentzVector(0, 0, 0, ParticleTypeDatabase::Proton.Mass());
        const TLorentzVector missing = beam_target - ggg;

        MM = missing.M();

        tree->Fill();
    }

}

double OmegaEtaG2::calcEnergySum2(const Event::Data &e) const
{
    double Sum = 0.0;

    if(mode != DataMode::MCTrue) {
        for(const auto& c : e.Candidates()) {
            const auto d = geo.DetectorFromAngles(c->Theta(),c->Phi());
            if( d & Detector_t::Any_t::CB) {
                Sum += c->ClusterEnergy();
            }

        }
    }

    for(const auto& c : e.Particles().GetAll()) {
        const auto d = geo.DetectorFromAngles(c->Theta(),c->Phi());
        if( d & Detector_t::Any_t::CB) {
            Sum += c->Ek();
        }

    }
    return Sum;
}

OmegaEtaG2::OmegaEtaG2(PhysOptPtr opts):
    OmegaBase("OmegaEtaG2", opts)
{
    if(opts->GetOption("Proton") == "MCTrue") {
        data_proton = false;
        LOG(INFO) << "Using proton from MCTrue";
    };

    if(opts->GetOption("Tagger") == "MCTrue") {
        data_tagger = false;
        LOG(INFO) << "Using Tagger from MCTrue";
    };

    if(opts->GetOption("ESum") != "") {
        ESum_cut = atof(opts->GetOption("ESum").c_str());
    }

    tree = new TTree("omegaetag2","");

    tree->Branch("pEk",     &pEk);
    tree->Branch("pTheta",  &pTheta);
    tree->Branch("pPhi",    &pPhi);
    tree->Branch("gggIM",   &gggIM);
    tree->Branch("gggTheta",&gggTheta);
    tree->Branch("gggPhi",  &gggPhi);
    tree->Branch("gggTime", &gggTime);
    tree->Branch("ggIM",     ggIM, "ggIM[3]/D");
    tree->Branch("MM",      &MM);
    tree->Branch("tagch",   &tagch);
    tree->Branch("tagtime", &tagtime);
    tree->Branch("rf",      &rf);
}

OmegaEtaG2::~OmegaEtaG2()
{

}

AUTO_REGISTER_PHYSICS(OmegaEtaG, "OmegaEtaG")
AUTO_REGISTER_PHYSICS(OmegaMCTruePlots, "OmegaMCTruePlots")
AUTO_REGISTER_PHYSICS(OmegeMCTree, "OmegaMCTree")
AUTO_REGISTER_PHYSICS(OmegaEtaG2, "OmegaEtaG2")
