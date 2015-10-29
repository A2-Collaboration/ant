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
#include "base/ParticleTypeTree.h"

#include "utils/particle_tools.h"
#include "utils/matcher.h"

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


OmegaEtaG::OmegaEtaG(const std::string& name, PhysOptPtr opts):
    OmegaBase(name, opts)
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

    auto pref = utils::ParticleTools::SanitizeDecayString(title);

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

    const string& decaystring = utils::ParticleTools::GetDecayString(event.MCTrue().ParticleTree());

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



OmegaMCTruePlots::OmegaMCTruePlots(const std::string& name, PhysOptPtr opts):
    Physics(name, opts)
{

}

void OmegaMCTruePlots::ProcessEvent(const Event& event)
{
    const auto& decaystring =utils::ParticleTools::GetProductionChannelString(event.MCTrue().ParticleTree());

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



TLorentzVector OmegaMCTree::getGamma1() const
{
    return gamma1_vector;
}

void OmegaMCTree::setGamma1(const TLorentzVector& value)
{
    gamma1_vector = value;
}

OmegaMCTree::OmegaMCTree(const std::string& name, PhysOptPtr opts): Physics(name, opts) {
    tree=new TTree("omegatree","omgega eta gamma MC true");
    tree->Branch("p", &proton_vector);
    tree->Branch("omega", &omega_vector);
    tree->Branch("gamma1", &gamma1_vector);
    tree->Branch("eta", &eta_vector);
    tree->Branch("gamma2", &gamma2_vector);
    tree->Branch("gamma3", &gamma3_vector);
}

OmegaMCTree::~OmegaMCTree()
{

}

void OmegaMCTree::ProcessEvent(const Event& event)
{
    if(!event.MCTrue().ParticleTree())
        return;

    struct TreeItem_t {
        const ParticleTypeDatabase::Type& Type;
        TLorentzVector* LorentzVector;
        TreeItem_t(const ParticleTypeDatabase::Type& type,
                   TLorentzVector* lv
                   ) :
            Type(type),
            LorentzVector(lv)
        {}
        // this operator makes Tree::Sort work
        bool operator<(const TreeItem_t& rhs) const {
            return Type.Name() < rhs.Type.Name();
        }
    };

    auto signal_tree = Tree<TreeItem_t>::MakeNode(ParticleTypeDatabase::BeamProton, (TLorentzVector*) nullptr);
    signal_tree->CreateDaughter(ParticleTypeDatabase::Proton, &proton_vector);
    auto omega = signal_tree->CreateDaughter(ParticleTypeDatabase::Omega, &omega_vector);
    omega->CreateDaughter(ParticleTypeDatabase::Photon, &gamma1_vector);
    auto eta = omega->CreateDaughter(ParticleTypeDatabase::Eta, &eta_vector);
    eta->CreateDaughter(ParticleTypeDatabase::Photon, &gamma2_vector);
    eta->CreateDaughter(ParticleTypeDatabase::Photon, &gamma3_vector);

    signal_tree->Sort();

    auto comparer = [] (const ParticlePtr& p, const TreeItem_t& item) {
        if(p->Type().Name() == item.Type.Name()) {
            if(item.LorentzVector)
                *item.LorentzVector = *p;
            return true;
        }
        return false;
    };

    if(event.MCTrue().ParticleTree()->IsEqual(signal_tree, comparer))
        tree->Fill();
}

void OmegaMCTree::ShowResult()
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
        v += p->Candidate()->Time();
        ++begin;
    }

    return v;
}

TLorentzVector boost(const TLorentzVector& lv, const TVector3& boot) {
    TLorentzVector v(lv);
    v.Boost(boot);
    return v;
}

struct chi2_highscore_t {

    double chi2  = std::numeric_limits<double>::infinity();
    int    index = -1;

    chi2_highscore_t() {};

    void Put(const double newchi2, const int newindex) {
        if(newchi2 < chi2) {
            chi2 = newchi2;
            index = newindex;
        }
    }
};

void OmegaEtaG2::Analyse(const Event::Data &data, const Event &event)
{

    const auto& particletree = event.MCTrue().ParticleTree();

    if(particletree) {

        if(utils::ParticleTools::FindParticle(ParticleTypeDatabase::Omega, particletree, 1)) {

            h_TotalEvents->Fill("#omega", 1);

            if(particletree->IsEqual(signal_tree, utils::ParticleTools::MatchByParticleName)) {
                h_TotalEvents->Fill("Signal",1);
            }
            else if(particletree->IsEqual(reference_tree, utils::ParticleTools::MatchByParticleName)) {
                h_TotalEvents->Fill("Reference",1);
            }
        }



    }

    const ParticleList& iphotons = data.Particles().Get(ParticleTypeDatabase::Photon);
    const ParticleList& iprotons = (data_proton ? event.Reconstructed() : event.MCTrue()).Particles().Get(ParticleTypeDatabase::Proton);

    steps->Fill("0 Events seen", 1);

    const auto Esum = data.TriggerInfos().CBEenergySum();

    if(Esum < 550.0)
        return;

    steps->Fill("1 CBEsum", 1);

    if(iphotons.size() != 3)
        return;

    if(iprotons.size() != 1)
        return;


    const ParticleList protons = FilterParticles(getGeoAccepted(iprotons), proton_cut);

    if(protons.size() != 1)
        return;

    const ParticleList photons = FilterParticles(getGeoAccepted(iphotons), photon_cut);

    if(photons.size() != 3)
        return;

    steps->Fill("2 nPhotons nProtons", 1);

    const auto& proton = protons.at(0);



    //const auto& mctrue_photons = event.MCTrue().Particles().Get(ParticleTypeDatabase::Photon);
    rf = static_cast<int>(identify(event));

    pbranch = mParticleVars(*proton);
    pTime  = data_proton ? proton->Candidate()->Time() : 0.0;

    const Particle ggg(ParticleTypeDatabase::Omega, LVSum(photons.begin(), photons.end()));
    gggTime  = TimeAvg(photons.begin(), photons.end());
    gggbranch = mParticleVars(ggg);


    const auto compl_angle = fabs(TVector2::Phi_mpi_pi(proton->Phi() - ggg.Phi()));

    if(!complcut.Contains(compl_angle))
        return;

    steps->Fill("3 Coplanarity", 1);

    g1branch = mParticleVars(*photons.at(0));
    g2branch = mParticleVars(*photons.at(1));
    g3branch = mParticleVars(*photons.at(2));

    const TVector3 gggBoost = -ggg.BoostVector();

    Chi2_Omega = std_ext::sqr((gggbranch.IM - omega_peak.Mean) / omega_peak.Sigma);




    for(const TaggerHitPtr& t : data_tagger?data.TaggerHits():event.MCTrue().TaggerHits()) {

        tagch   = t->Channel();
        tagtime = t->Time();

        const TLorentzVector beam_target = t->PhotonBeam() + TLorentzVector(0, 0, 0, ParticleTypeDatabase::Proton.Mass());
        const Particle missing(ParticleTypeDatabase::Proton, beam_target - ggg);

        calcp = analysis::utils::ParticleVars(missing);

        angle_p_calcp = radian_to_degree(missing.Angle(proton->Vect()));

        int combindex = 0;

        chi2_highscore_t best_eta;
        chi2_highscore_t best_pi0;

        const vector<vector<size_t>> combs = {{0,1,2},{0,2,1},{1,2,0}};
        for(const auto& comb : combs) {
            const auto& g1 = photons.at(comb[0]);
            const auto& g2 = photons.at(comb[1]);
            const auto& g3 = photons.at(comb[2]);

            const TLorentzVector gg = *g1 + *g2;

            ggIM[combindex] = gg.M();

            Chi2_Pi0[combindex] =  std_ext::sqr((ggIM[combindex] - pi0_peak.Mean) / pi0_peak.Sigma);
            Chi2_Eta[combindex] =  std_ext::sqr((ggIM[combindex] - eta_peak.Mean) / eta_peak.Sigma);

            const TLorentzVector g3_boosted = boost(*g3, gggBoost);

            EgOmegaSys[combindex] = g3_boosted.E();

            best_eta.Put(Chi2_Eta[combindex], combindex);
            best_pi0.Put(Chi2_Pi0[combindex], combindex);

            ++combindex;

        }

        bestEtaIn = best_eta.index;
        bestPi0In = best_pi0.index;

        if(best_eta.chi2 < best_pi0.chi2) {
            bestChi = best_eta.chi2;
            bestHyp = 1;
        } else {
            bestChi = best_pi0.chi2;
            bestHyp = 2;
        }

        tree->Fill();

    }

}

OmegaEtaG2::channel_type_t OmegaEtaG2::identify(const Event &event) const
{

    auto particletree = event.MCTrue().ParticleTree();

    if(!particletree)
        return channel_type_t::Background;

    if(particletree->IsEqual(signal_tree, utils::ParticleTools::MatchByParticleName))
        return channel_type_t::Signal;
    else if(particletree->IsEqual(reference_tree, utils::ParticleTools::MatchByParticleName))
        return channel_type_t::Reference;
    else
        return channel_type_t::Background;
}

ParticleList OmegaEtaG2::FilterParticles(const data::ParticleList& list, const particleCuts_t& cuts) const {
    ParticleList olist;
    //copy_if(list.begin(), list.end(), olist.begin(), [cuts] (const ParticlePtr& p) { return cuts.TestParticle(*p);});
    for(const auto& p : list) {
        if(cuts.TestParticle(*p)) {
            olist.emplace_back(p);
        }
    }
    return olist;
}

OmegaEtaG2::OmegaEtaG2(const std::string& name, PhysOptPtr opts):
    OmegaBase(name, opts)
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

    proton_cut.E_range     = { 50,1000 };
    proton_cut.Theta_range = degree_to_radian(interval<double>( 2, 45));

    photon_cut.E_range     = { 0, 1600 };
    photon_cut.Theta_range = degree_to_radian(interval<double>( 2, 160));


    tree = HistFac.makeTTree("tree");

    pbranch.SetBranches(tree, "p");
    tree->Branch("pTime", &pTime);

    gggbranch.SetBranches(tree, "ggg");
    tree->Branch("gggTime", &gggTime);

    tree->Branch("AnglePcP", &angle_p_calcp);

    tree->Branch("ggIM",     ggIM, "ggIM[3]/D");

    calcp.SetBranches(tree,"calcp");

    g1branch.SetBranches(tree, "g1");
    g2branch.SetBranches(tree, "g2");
    g3branch.SetBranches(tree, "g3");

    tree->Branch("tagch",   &tagch);
    tree->Branch("tagtime", &tagtime);
    tree->Branch("rf",      &rf);

    tree->Branch("chi2_omega",    &Chi2_Omega);
    tree->Branch("chi2_eta[3]",      Chi2_Eta,"chi2_eta[3]/D");
    tree->Branch("chi2_pi0[3]",      Chi2_Pi0,"chi2_pi0[3]/D");
    tree->Branch("EgOmegaSys[3]",    EgOmegaSys,"EgOmegaSys[3]/D");
    tree->Branch("ibestEta", &bestEtaIn);
    tree->Branch("ibestPi0", &bestPi0In);
    tree->Branch("bestChi",  &bestChi);
    tree->Branch("fbestHyp",  &bestHyp);

    signal_tree = ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Omega_gEta_3g);
    reference_tree = ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Omega_gPi0_3g);

    steps = HistFac.makeTH1D("Steps","Step","Events passed",BinSettings(14),"steps");

    h_TotalEvents = HistFac.makeTH1D("TotalEvents","","",BinSettings(3),"TotalEvents");


}

OmegaEtaG2::~OmegaEtaG2()
{

}

void OmegaEtaG2::mParticleVars::SetBranches(TTree* tree, const string& name)
{
    tree->Branch((name+"matchAngle").c_str(), &matchAngle);
    ParticleVars::SetBranches(tree,name);
}



bool OmegaEtaG2::particleCuts_t::TestParticle(const Particle& p) const
{
    return E_range.Contains(p.Ek()) && Theta_range.Contains(p.Theta());
}

AUTO_REGISTER_PHYSICS(OmegaEtaG)
AUTO_REGISTER_PHYSICS(OmegaMCTruePlots)
AUTO_REGISTER_PHYSICS(OmegaMCTree)
AUTO_REGISTER_PHYSICS(OmegaEtaG2)
