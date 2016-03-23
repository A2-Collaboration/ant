#include "omega.h"
#include "TH1D.h"
#include "plot/root_draw.h"
#include "utils/combinatorics.h"
#include <string>
#include <iostream>
#include "TH3.h"
#include "base/Logger.h"
#include <algorithm>
#include <iostream>
#include "base/std_ext/math.h"

#include "TTree.h"
#include "base/std_ext/iterators.h"
#include "base/ParticleTypeTree.h"

#include "utils/particle_tools.h"
#include "utils/matcher.h"

#include "APLCON.hpp"
#include "expconfig/ExpConfig.h"
#include "base/WrapTFile.h"
#include "TCanvas.h"
#include <cassert>

#include "utils/matcher.h"

#include "root-addons/analysis_codes/hstack.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace ant::std_ext;

void OmegaBase::ProcessEvent(const TEvent& event, manager_t& manager)
{
    const auto& data = mode==DataMode::Reconstructed ? event.Reconstructed() : event.MCTrue();
    Analyse(data, event, manager);
}

double OmegaBase::calcEnergySum(const TParticleList &particles) const
{
    double esum = 0.0;

    for( const TParticlePtr& p : particles) {
        if( geo.DetectorFromAngles(p->Theta(), p->Phi()) == Detector_t::Type_t::CB ) {
            esum += p->Ek();
        }
    }

    return esum;
}

TParticleList OmegaBase::getGeoAccepted(const TParticleList &p) const
{
    TParticleList list;
    for( auto& particle : p) {
        if( geo.DetectorFromAngles(particle->Theta(), particle->Phi()) != Detector_t::Any_t::None )
            list.emplace_back(particle);
    }
    return list;
}

unsigned OmegaBase::geoAccepted(const TCandidateList& cands) const {

    unsigned n = 0;

    for( auto& c : cands) {
        if( geo.DetectorFromAngles(c.Theta, c.Phi) != Detector_t::Any_t::None )
            ++n;
    }

    return n;
}

OmegaBase::OmegaBase(const string &name, OptionsPtr opts):
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


OmegaEtaG::OmegaEtaG(const std::string& name, OptionsPtr opts):
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

    if(Options->Get<string>("OmegaEtaGMode") == "McTrue") {
        mode = OmegaBase::DataMode::MCTrue;
    }

    VLOG(8) << "mode option is " << Options->Get<string>("OmegaEtaGMode");
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

LorentzVec Boost(const LorentzVec& lv, const vec3& boost) {
    LorentzVec b(lv);
    b.Boost(boost);
    return b;
}

void OmegaEtaG::Analyse(const TEventData &data, const TEvent& event, manager_t&)
{

    steps->Fill("Events seen",1);

    const auto nPhotons = data.Particles.Get(ParticleTypeDatabase::Photon).size();
    const auto nProtons = data.Particles.Get(ParticleTypeDatabase::Proton).size();

    perDecayhists_t* h = nullptr;

    const string& decaystring = utils::ParticleTools::GetDecayString(event.MCTrue().ParticleTree);

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

    if( data.Candidates.size() > 4 || data.Candidates.size() < 3 ) // not more then that (3photns +0or1 protons)
        return;

    steps->Fill("nCandidates",1);

    const double CBESum = mode==DataMode::Reconstructed ? calcEnergySum(data.Particles.Get(ParticleTypeDatabase::Photon)) : data.Trigger.CBEnergySum;

    if( CBESum < 550.0 )
        return;

    steps->Fill("ESum",1);

    const TParticleList photons = getGeoAccepted(data.Particles.Get(ParticleTypeDatabase::Photon));
    const TParticleList protons = getGeoAccepted(data.Particles.Get(ParticleTypeDatabase::Proton));


    bool is_omega_decay = false;
    const ParticleTypeDatabase::Type* subtype = nullptr;



    for( auto comb = utils::makeCombination(photons,3); !comb.Done(); ++comb) {

        TParticleList ggg_list;
        ggg_list.assign(comb.begin(),comb.end());

        LorentzVec gggState = *comb.at(0)+*comb.at(1)+*comb.at(2);
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

        const auto& mc_protons = event.MCTrue().Particles.Get(ParticleTypeDatabase::Proton);

        const TParticlePtr mc_p = !mc_protons.empty() ? mc_protons.at(0) : shared_ptr<TParticle>(nullptr);

        if(omega_range.Contains(gggIM) && h) {
            for(auto& th : event.MCTrue().TaggerHits) {
                const auto beam_target = th.GetPhotonBeam() + LorentzVec(0, 0, 0, ParticleTypeDatabase::Proton.Mass());
                const auto mm = beam_target - gggState;
                const auto mm_boosted = Boost(mm, -beam_target.BoostVector());


                h->calc_proton_energy_theta->Fill(mm.E-ParticleTypeDatabase::Proton.Mass(), mm.Theta() * TMath::RadToDeg());

                h->mm->Fill(mm.M());

                if(mc_p) {

                    const TParticle p_special(ParticleTypeDatabase::Proton, mm.E-ParticleTypeDatabase::Proton.Mass(), mc_p->Theta(), mc_p->Phi());
                    h->calc_proton_special->Fill(p_special.Ek(), p_special.Theta()*TMath::RadToDeg());

                    LorentzVec mc_p_v = *mc_p;
                    mc_p_v.Boost(-beam_target.BoostVector());
                    const auto angle = mc_p_v.Angle(mm_boosted.p);
                    h->angle_p->Fill(angle * TMath::RadToDeg());

                    LorentzVec ggg_boost = gggState;
                    ggg_boost.Boost(-beam_target.BoostVector());
                    const auto angle_pggg = ggg_boost.Angle(mc_p_v.p);
                    h->angle_p_ggg->Fill(angle_pggg * TMath::RadToDeg());

                    const auto p_phi_diff = vec2::Phi_mpi_pi(mm_boosted.Phi() - mc_p_v.Phi());
                    h->p_phi_diff->Fill(p_phi_diff);

                }
            }
        }

        for( auto gcomb = utils::makeCombination(ggg_list,2); !gcomb.Done(); ++gcomb) {

            const LorentzVec g1(*gcomb.at(0));
            const LorentzVec g2(*gcomb.at(1));
            const LorentzVec ggState = g1 + g2;
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


    auto cit = getCircularIterator(ColorPalette::Colors.begin(), ColorPalette::Colors.end());

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

    canvas("OmegaEtaG per Decay Results gg") << &stack << endc;
    canvas("OmegaEtaG per Decay Results ggg") << &stack2 << endc;
    canvas("OmegaEtaG per Decay Results mm") << drawoption("pads") << &stack3 << endc;
    canvas("OmegaEtaG per Decay Results Angle p") << drawoption("pads") << &stack4 << endc;
    canvas("OmegaEtaG per Decay Results Angle p ggg") << drawoption("pads") << &stack5 << endc;
    canvas("OmegaEtaG per Decay Resultsp phi diff") << drawoption("pads") << &stack6 << endc;
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




OmegaMCTruePlots::PerChannel_t::PerChannel_t(const string& Title, HistogramFactory& hf):
    title(Title)
{
    proton_E_theta = hf.makeTH2D(title,"E [MeV]","#theta [#circ]",BinSettings(1000),BinSettings(360,0,180), title+"_e_theta");
}

void OmegaMCTruePlots::PerChannel_t::Show()
{
    canvas("Omega per Channel: "+title) << drawoption("colz") << proton_E_theta << endc;
}

void OmegaMCTruePlots::PerChannel_t::Fill(const TEventData& d)
{
    const auto& protons = d.Particles.Get(ParticleTypeDatabase::Proton);
    if(!protons.empty()) {
        const auto& p = protons.at(0);
        proton_E_theta->Fill(p->Ek(), p->Theta()*TMath::RadToDeg());
    }
}



OmegaMCTruePlots::OmegaMCTruePlots(const std::string& name, OptionsPtr opts):
    Physics(name, opts)
{

}

void OmegaMCTruePlots::ProcessEvent(const TEvent& event, manager_t&)
{
    const auto& decaystring =utils::ParticleTools::GetProductionChannelString(event.MCTrue().ParticleTree);

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



LorentzVec OmegaMCTree::getGamma1() const
{
    return gamma1_vector;
}

void OmegaMCTree::setGamma1(const LorentzVec& value)
{
    gamma1_vector = value;
}

OmegaMCTree::OmegaMCTree(const std::string& name, OptionsPtr opts): Physics(name, opts) {
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

void OmegaMCTree::ProcessEvent(const TEvent& event, manager_t&)
{
    if(!event.MCTrue().ParticleTree)
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

    auto comparer = [] (const TParticlePtr& p, const TreeItem_t& item) {
        if(p->Type().Name() == item.Type.Name()) {
            if(item.LorentzVector)
                *item.LorentzVector = *p;
            return true;
        }
        return false;
    };

    if(event.MCTrue().ParticleTree->IsEqual(signal_tree, comparer))
        tree->Fill();
}

void OmegaMCTree::ShowResult()
{

}


template <typename it_type>
LorentzVec LVSum(it_type begin, it_type end) {
    LorentzVec v;

    while(begin!=end) {
        v += **begin;
        ++begin;
    }

    return v;
}

template <typename it_type>
LorentzVec LVSumL(it_type begin, it_type end) {
    LorentzVec v;

    while(begin!=end) {
        v += *begin;
        ++begin;
    }

    return v;
}

template <typename it_type>
double TimeAvg(it_type begin, it_type end) {
    double t    = 0.0;
    double Esum = 0.0;

    while(begin!=end) {
        const TParticlePtr& p = *begin;
        t += p->Candidate->Time * p->Candidate->CaloEnergy;
        Esum += p->Candidate->CaloEnergy;
        ++begin;
    }

    return t / Esum;
}


#define FASSERT(x) if(!(x)) LOG(ERROR) << "ERROR";

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


double IM(const TParticlePtr& p1, const TParticlePtr& p2) {
    return (*p1+*p2).M();
}

double getTime(const TParticlePtr& p) {
    return p->Candidate != nullptr ? p->Candidate->Time : std_ext::NaN;
}

void OmegaEtaG2::Analyse(const TEventData &data, const TEvent& event, manager_t& manager)
{

    const auto& particletree = event.MCTrue().ParticleTree;

    //const auto& mctrue_photons = event.MCTrue().Particles().Get(ParticleTypeDatabase::Photon);
    t.Channel = reaction_channels.identify(event.MCTrue().ParticleTree);

    if(t.Channel == ReactionChannelList_t::other_index) {
        if(event.MCTrue().ParticleTree!=nullptr) {
            missed_channels->Fill(utils::ParticleTools::GetDecayString(event.MCTrue().ParticleTree).c_str(), 1.0);
        }
    } else {
        found_channels->Fill(t.Channel);
    }

    TH1* steps = stephists.at(t.Channel);

    steps->Fill("0 Events seen", 1);

    const auto Esum = data.Trigger.CBEnergySum;

    if(Esum <  cut_ESum)
        return;

    steps->Fill("1 CBEsum", 1);

    const auto n_cands = geoAccepted(data.Candidates);

    if(n_cands != 4) {
        return;
    }

    steps->Fill("2 nCands", 1);


    TParticleList iphotons;
    TParticleList iprotons;

    for(auto p: data.Candidates.get_iter()) {
        if(p->VetoEnergy < .25) {
            iphotons.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::Photon, p));
        } else {
            iprotons.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::Proton, p));
        }
    }

    if(iphotons.size() != 3)
        return;

    if(iprotons.size() != 1)
        return;

    const TParticleList protons = FilterProtons(getGeoAccepted(iprotons));

    if(protons.size() != 1)
        return;

    const TParticleList photons = FilterPhotons(getGeoAccepted(iphotons));

    if(photons.size() != 3)
        return;

    steps->Fill("3 nPhotons nProtons", 1);

    const auto& proton = protons.at(0);

    t.photons().at(0) = *photons.at(0);
    t.photons().at(1) = *photons.at(1);
    t.photons().at(2) = *photons.at(2);

    t.p      = *proton;
    t.p_Time = getTime(proton);

    t.p_PSA_Angle  = std_ext::NaN;
    t.p_PSA_Radius = std_ext::NaN;
    t.p_detector   = 0;

    if(proton->Candidate) {

        if(proton->Candidate->Detector & Detector_t::Type_t::TAPS) {
            t.p_detector = 2;
            const auto& cluster = proton->Candidate->FindCaloCluster();
            if(cluster) {
                t.p_PSA_Angle  = radian_to_degree(cluster->GetPSAAngle());
                t.p_PSA_Radius = cluster->GetPSARadius();
            }
        } else if(proton->Candidate->Detector & Detector_t::Type_t::CB) {
            t.p_detector = 1;
        }

    }


    const TParticle ggg(ParticleTypeDatabase::Omega, LVSum(photons.begin(), photons.end()));
    t.ggg = ggg;
    const auto gggBoost = -ggg.BoostVector();

    t.copl_angle = fabs(vec2::Phi_mpi_pi(proton->Phi() - ggg.Phi() - M_PI));

    if(t.copl_angle > cut_Copl)
        return;

    steps->Fill("4 Coplanarity", 1);

    t.CBAvgTime = event.Reconstructed().Trigger.CBTiming;
    if(!isfinite(t.CBAvgTime))
        return;

    steps->Fill("5 valid CB Avg Time", 1);

    if(data.TaggerHits.size() > 0)
        steps->Fill("6 has TaggHits", 1);

    for(const TTaggerHit& TagH : data.TaggerHits) {

        promptrandom.SetTaggerHit(TagH.Time - t.CBAvgTime);

        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        t.TaggW  = promptrandom.FillWeight();
        t.TaggE  = TagH.PhotonEnergy;
        t.TaggCh = TagH.Channel;
        t.TaggT  = TagH.Time;

        const LorentzVec beam_target = TagH.GetPhotonBeam() + LorentzVec(0, 0, 0, ParticleTypeDatabase::Proton.Mass()); // make global
        const TParticle missing(ParticleTypeDatabase::Proton, beam_target - ggg);

        t.mm = missing;

        t.p_mm_angle = radian_to_degree(missing.Angle(*proton));


        fitter.SetEgammaBeam(TagH.PhotonEnergy);
        fitter.SetProton(proton);
        fitter.SetPhotons(photons);

        auto fitres = fitter.DoFit();

        if(fitres.Status != APLCON::Result_Status_t::Success)
            continue;

        t.KinFitChi2 = fitres.ChiSquare / fitres.NDoF;
        t.KinFitIterations = unsigned(fitres.NIterations);

        t.p_fitted = *fitter.GetFittedProton();

        t.photons_fitted().at(0) = *fitter.GetFittedPhotons().at(0);
        t.photons_fitted().at(1) = *fitter.GetFittedPhotons().at(1);
        t.photons_fitted().at(2) = *fitter.GetFittedPhotons().at(2);

        const TParticle ggg_fitted(ParticleTypeDatabase::Omega, LVSumL(t.photons_fitted().begin(), t.photons_fitted().end()));
        t.ggg_fitted = ggg_fitted;

        const auto gggBoost_fitted = -ggg_fitted.BoostVector();

        size_t combindex = 0;

        for(const auto& comb : combs) {
            const auto& g1 = photons.at(comb[0]);
            const auto& g2 = photons.at(comb[1]);
            const auto& g3 = photons.at(comb[2]);

            const auto gg = *g1 + *g2;

            t.ggIM().at(combindex) = gg.M();

            const auto g3_boosted = Boost(*g3, gggBoost);

            t.BachelorE().at(combindex) = g3_boosted.E;

            ++combindex;
        }

        combindex = 0;
        for(const auto& comb : combs) {
            const auto& g1 = t.photons_fitted().at(comb[0]);
            const auto& g2 = t.photons_fitted().at(comb[1]);
            const auto& g3 = t.photons_fitted().at(comb[2]);

            const auto gg = g1 + g2;

            t.ggIM_fitted().at(combindex) = gg.M();

            const auto g3_boosted = Boost(g3, gggBoost_fitted);

            t.BachelorE_fitted().at(combindex) = g3_boosted.E;

            ++combindex;
        }




        TParticleList rec_photons(3);
        TParticlePtr  rec_proton = nullptr;
        TParticleList true_particles(4);

        t.ggIM_real = NaN;
        t.ggIM_comb = {NaN, NaN};

        if(particletree && (t.Channel == 1 || t.Channel == 2)) {
            particletree->Map_level([&true_particles] (const TParticlePtr& p, const size_t& level) {

                if(level == 1) {
                    if(p->Type() == ParticleTypeDatabase::Proton) {
                        FASSERT(true_particles[3] == nullptr);
                        true_particles[3] = p;
                    }
                }

                if(p->Type() == ParticleTypeDatabase::Photon) {

                    if(level==2) {
                        FASSERT(true_particles[0]==nullptr);
                        true_particles[0] = p;

                    } else if(level == 3) {

                        if(!true_particles[1]) {
                            true_particles[1] = p;
                        } else {
                            FASSERT(true_particles[2]==nullptr);
                            true_particles[2] = p;
                        }
                    }
                }
            });

            FASSERT(true_particles[0]!=nullptr);
            FASSERT(true_particles[1]!=nullptr);
            FASSERT(true_particles[2]!=nullptr);
            FASSERT(true_particles[3]!=nullptr);

            t.p_true = *true_particles[3];

            const auto matched  = utils::match1to1(true_particles, data.Particles.GetAll(), TParticle::CalcAngle, {0.0, degree_to_radian(15.0)});

            if(matched.size() == true_particles.size()) {

                rec_photons[0] = utils::FindMatched(matched, true_particles[0]);
                rec_photons[1] = utils::FindMatched(matched, true_particles[1]);
                rec_photons[2] = utils::FindMatched(matched, true_particles[2]);
                rec_proton     = utils::FindMatched(matched, true_particles[3]);

                FASSERT(rec_photons[0]!=nullptr);
                FASSERT(rec_photons[1]!=nullptr);
                FASSERT(rec_photons[2]!=nullptr);
                FASSERT(rec_proton    !=nullptr);


                t.p_matched = (rec_proton->Type() == ParticleTypeDatabase::Proton);

                t.ggIM_real      = IM(rec_photons[1], rec_photons[2]);
                t.ggIM_comb()[0] = IM(rec_photons[0], rec_photons[1]);
                t.ggIM_comb()[1] = IM(rec_photons[0], rec_photons[2]);

            }

        }

        if(opt_save_after_kinfit && t.KinFitChi2 < 5.0)
            manager.SaveEvent();

        tree->Fill();

    }

}

OmegaEtaG2::ReactionChannelList_t OmegaEtaG2::makeChannels()
{
    ReactionChannelList_t m;

    m.channels[0] = ReactionChannel_t("Data");
    m.channels[1] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Omega_gEta_3g), kRed};  //sig
    m.channels[2] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Omega_gPi0_3g), kGreen};  //ref
    m.channels[3] = ReactionChannel_t("Sum MC");
    m.channels[4] = ReactionChannel_t("MC BackG");

    m.channels[10] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0_2g),      "#pi^{0}", kYellow};
    m.channels[11] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::TwoPi0_4g),   "#pi^{0} #pi^{0}", kOrange};
    m.channels[12] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::ThreePi0_6g), "#pi^{0} #pi^{0} #pi^{0}",kGreen-9};
    m.channels[13] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0Eta_4g),   "#pi^{0} #eta",kBlue};
    m.channels[14] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Omega_Pi0PiPPiM_2g),"#omega #rightarrow #pi^{0} #pi^{+} #pi^{-}",kMagenta};
    m.channels[m.other_index] = ReactionChannel_t(nullptr, "Others", kCyan);

    return m;
}

bool OmegaEtaG2::AcceptedPhoton(const TParticlePtr& photon)
{
    if(photon->Candidate->Detector & Detector_t::Type_t::CB) {
        if(photon_E_cb.Contains(photon->Ek())) {
            return true;
        }
    } else if(photon->Candidate->Detector & Detector_t::Type_t::TAPS) {
        if(photon_E_taps.Contains(photon->Ek())) {
            return true;
        }
    }

    return false;
}

bool OmegaEtaG2::AcceptedProton(const TParticlePtr& proton)
{
    if(proton_theta.Contains(proton->Theta())){
        return true;
    }

    return false;
}

TParticleList OmegaEtaG2::FilterPhotons(const TParticleList& list)
{

    TParticleList olist;

    for(const auto& p : list) {
        if(AcceptedPhoton(p)) {
            olist.emplace_back(p);
        }
    }
    return olist;
}

TParticleList OmegaEtaG2::FilterProtons(const TParticleList& list)
{

    TParticleList olist;

    for(const auto& p : list) {
        if(AcceptedProton(p)) {
            olist.emplace_back(p);
        }
    }
    return olist;
}


OmegaEtaG2::OmegaEtaG2(const std::string& name, OptionsPtr opts):
    OmegaBase(name, opts),
    tree(HistFac.makeTTree("tree")),

    cut_ESum(opts->Get<double>("CBESum", 550.0)),
    cut_Copl(degree_to_radian(opts->Get<double>("CoplAngle", 15.0))),
    photon_E_cb(opts->Get<decltype(photon_E_cb)>("PhotonECB", {50.0, 1600.0})),
    photon_E_taps(opts->Get<decltype(photon_E_taps)>("PhotonETAPS", {200.0, 1600.0})),
    proton_theta(degree_to_radian(opts->Get<decltype(proton_theta)>("ProtonThetaRange", {2.0, 45.0}))),
    fitter("OmegaEtaG2", 3, utils::UncertaintyModels::MCExtracted::makeAndLoad())
{

    promptrandom.AddPromptRange({-5,5});
    promptrandom.AddRandomRange({-20, -10});
    promptrandom.AddRandomRange({ 10,  20});

    t.CreateBranches(tree);

    missed_channels = HistFac.makeTH1D("Unlisted Channels","","Total Events seen",BinSettings(20),"unlistedChannels");
    found_channels  = HistFac.makeTH1D("Listed Channels",  "","Total Events seen",BinSettings(20),"listedChannels");

    for(const auto& c : reaction_channels.channels) {

        stephists[c.first] = HistFac.makeTH1D("Steps: " + c.second.name, "", "", BinSettings(14), "steps_" + to_string(c.first));
        stephists[c.first]->SetLineColor(c.second.color);

        if(c.first<20)
            found_channels->GetXaxis()->SetBinLabel(c.first+1,c.second.name.c_str());
    }

    opt_save_after_kinfit = opts->Get("SaveAfterKinfit", opt_save_after_kinfit);

}

OmegaEtaG2::~OmegaEtaG2()
{
}

void OmegaEtaG2::Finish()
{
    hstack* s = HistFac.make<hstack>("steps");
    for(auto& c : stephists) {
        (*s) << c.second;
    }
}



OmegaEtaG2::OmegaTree_t::OmegaTree_t()
{}

decltype(OmegaEtaG2::combs) OmegaEtaG2::combs = {{0,1,2},{0,2,1},{1,2,0}};



OmegaEtaG2::ReactionChannel_t::ReactionChannel_t(const std::shared_ptr<OmegaEtaG2::decaytree_t> &t, const string &n, const int c):
    name(n), tree(t), color(c)
{
}

OmegaEtaG2::ReactionChannel_t::ReactionChannel_t(const std::shared_ptr<OmegaEtaG2::decaytree_t> &t, const int c):
    name(utils::ParticleTools::GetDecayString(t)),
    tree(t),
    color(c)
{}

OmegaEtaG2::ReactionChannel_t::ReactionChannel_t(const string &n):
    name(n)
{}

OmegaEtaG2::ReactionChannel_t::~ReactionChannel_t()
{}



unsigned OmegaEtaG2::ReactionChannelList_t::identify(const ant::TParticleTree_t& tree) const
{

    if(!tree)
        return 0;

    for(const auto& c : channels) {

        if(!c.second.tree)
            continue;

        if(tree->IsEqual(c.second.tree, utils::ParticleTools::MatchByParticleName)) {
            return c.first;
        }
    }

    return other_index;
}

const OmegaEtaG2::ReactionChannelList_t OmegaEtaG2::reaction_channels = OmegaEtaG2::makeChannels();

const unsigned OmegaEtaG2::ReactionChannelList_t::other_index = 1000;

AUTO_REGISTER_PHYSICS(OmegaEtaG)
AUTO_REGISTER_PHYSICS(OmegaMCTruePlots)
AUTO_REGISTER_PHYSICS(OmegaMCTree)
AUTO_REGISTER_PHYSICS(OmegaEtaG2)
