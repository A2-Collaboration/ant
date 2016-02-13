#include "omega.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include "plot/root_draw.h"
#include "plot/Histogram.h"
#include "utils/combinatorics.h"
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

#include "APLCON.hpp"
#include "expconfig/ExpConfig.h"
#include "base/WrapTFile.h"
#include "TCanvas.h"
#include <cassert>

#include "utils/matcher.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace ant::std_ext;

void OmegaBase::ProcessEvent(const TEvent& event, manager_t& manager)
{
    const auto& data = mode==DataMode::Reconstructed ? *event.Reconstructed : *event.MCTrue;
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
        if( geo.DetectorFromAngles(c->Theta, c->Phi) != Detector_t::Any_t::None )
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

TLorentzVector Boost(const TLorentzVector& lv, const TVector3& boost) {
    TLorentzVector b(lv);
    b.Boost(boost);
    return b;
}

void OmegaEtaG::Analyse(const TEventData &data, const TEvent& event, manager_t&)
{

    steps->Fill("Events seen",1);

    const auto nPhotons = data.Particles.Get(ParticleTypeDatabase::Photon).size();
    const auto nProtons = data.Particles.Get(ParticleTypeDatabase::Proton).size();

    perDecayhists_t* h = nullptr;

    const string& decaystring = utils::ParticleTools::GetDecayString(event.MCTrue->ParticleTree);

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

        const auto& mc_protons = event.MCTrue->Particles.Get(ParticleTypeDatabase::Proton);

        const TParticlePtr mc_p = !mc_protons.empty() ? mc_protons.at(0) : shared_ptr<TParticle>(nullptr);

        if(omega_range.Contains(gggIM) && h) {
            for(auto& th : event.MCTrue->TaggerHits) {
                const TLorentzVector beam_target = th.GetPhotonBeam() + TLorentzVector(0, 0, 0, ParticleTypeDatabase::Proton.Mass());
                const TLorentzVector mm = beam_target - gggState;
                const TLorentzVector mm_boosted = Boost(mm,-beam_target.BoostVector());


                h->calc_proton_energy_theta->Fill(mm.Energy()-ParticleTypeDatabase::Proton.Mass(), mm.Theta() * TMath::RadToDeg());

                h->mm->Fill(mm.M());

                if(mc_p) {

                    const TParticle p_special(ParticleTypeDatabase::Proton, mm.Energy()-ParticleTypeDatabase::Proton.Mass(), mc_p->Theta(), mc_p->Phi());
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
    const auto& decaystring =utils::ParticleTools::GetProductionChannelString(event.MCTrue->ParticleTree);

    auto e = channels.find(decaystring);

    if(e == channels.end()) {
        channels.insert({decaystring, PerChannel_t(decaystring,HistFac)});
    }

    e = channels.find(decaystring);
    e->second.Fill(*event.MCTrue);
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
    if(!event.MCTrue->ParticleTree)
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

    if(event.MCTrue->ParticleTree->IsEqual(signal_tree, comparer))
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

TLorentzVector boost(const TLorentzVector& lv, const TVector3& boot) {
    TLorentzVector v(lv);
    v.Boost(boot);
    return v;
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
    return (TLorentzVector(*p1)+TLorentzVector(*p2)).M();
}


void OmegaEtaG2::Analyse(const TEventData &data, const TEvent& event, manager_t& manager)
{

    const auto& particletree = event.MCTrue->ParticleTree;

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

    if(just_preselect) {
        manager.SaveEvent();
        return;
    }

    TParticleList iphotons;
    TParticleList iprotons;

    for(auto p: data.Candidates) {
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


    //const auto& mctrue_photons = event.MCTrue().Particles().Get(ParticleTypeDatabase::Photon);
    branches.b_SigBgFlag = identify(event);

    branches.b_p = analysis::utils::ParticleVars(*proton);
    branches.b_pTime       = data_proton ? proton->Candidate->Time : numeric_limits<double>::quiet_NaN();

    branches.b_p_PSA_Angle = numeric_limits<double>::quiet_NaN();
    branches.b_p_PSA_R     = numeric_limits<double>::quiet_NaN();
    branches.b_p_detector  = 0;

    if(proton->Candidate) {

        if(proton->Candidate->Detector & Detector_t::Type_t::TAPS) {
            branches.b_p_detector = 2;
            const auto& cluster = proton->Candidate->FindCaloCluster();
            if(cluster) {
                branches.b_p_PSA_Angle = std_ext::radian_to_degree(cluster->GetPSAAngle());
                branches.b_p_PSA_R     = cluster->GetPSARadius();
            }
        } else if(proton->Candidate->Detector & Detector_t::Type_t::CB) {
            branches.b_p_detector = 1;
        }

    }


    const TParticle ggg(ParticleTypeDatabase::Omega, LVSum(photons.begin(), photons.end()));
    branches.b_gggTime  = TimeAvg(photons.begin(), photons.end());
    branches.b_ggg      = analysis::utils::ParticleVars(ggg);


    branches.b_copl_angle = fabs(TVector2::Phi_mpi_pi(proton->Phi() - ggg.Phi() - M_PI));

    if(branches.b_copl_angle > cut_Copl)
        return;

    steps->Fill("4 Coplanarity", 1);

    branches.b_g1 = analysis::utils::ParticleVars(*photons.at(0));
    branches.b_g2 = analysis::utils::ParticleVars(*photons.at(1));
    branches.b_g3 = analysis::utils::ParticleVars(*photons.at(2));

    const TVector3 gggBoost = -ggg.BoostVector();

    branches.b_CBAvgTime = event.Reconstructed->Trigger.CBTiming;
    if(!isfinite(branches.b_CBAvgTime))
        return;

    for(const TTaggerHit& t : data_tagger ? data.TaggerHits : event.MCTrue->TaggerHits) {

        promptrandom.SetTaggerHit(t.Time - branches.b_CBAvgTime);

        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        branches.b_TagW    = promptrandom.FillWeight();
        branches.b_TagE    = t.PhotonEnergy;
        branches.b_TagCh   = unsigned(t.Channel);
        branches.b_TagTime = t.Time;

        const TLorentzVector beam_target = t.GetPhotonBeam() + TLorentzVector(0, 0, 0, ParticleTypeDatabase::Proton.Mass());
        const TParticle missing(ParticleTypeDatabase::Proton, beam_target - ggg);

        branches.b_mmvector = analysis::utils::ParticleVars(missing);

        branches.b_p_mm_angle = radian_to_degree(missing.Angle(proton->Vect()));

        int combindex = 0;

        const vector<vector<size_t>> combs = {{0,1,2},{0,2,1},{1,2,0}};

        for(const auto& comb : combs) {
            const auto& g1 = photons.at(comb[0]);
            const auto& g2 = photons.at(comb[1]);
            const auto& g3 = photons.at(comb[2]);

            const TLorentzVector gg = *g1 + *g2;

            branches.b_ggIM[combindex] = gg.M();

            const TLorentzVector g3_boosted = boost(*g3, gggBoost);

            branches.b_BachelorE[combindex] = g3_boosted.E();

            ++combindex;

        }

        fitter.SetEgammaBeam(t.PhotonEnergy);
        fitter.SetProton(proton);
        fitter.SetPhotons(photons);

        auto fitres = fitter.DoFit();

        if(fitres.Status == APLCON::Result_Status_t::Success) {

            for(const auto& it_map : fitres.Variables) {
                const string& varname = it_map.first;
                const APLCON::Result_Variable_t& var = it_map.second;
                auto it_pull = pulls.find(varname);
                TH1D* h_pull;
                if(it_pull == pulls.end()) {
                    // not found so far, create the histogram on-the-fly
                    stringstream title;
                    title << "Pull " << var.PristineName << " " << component.at(var.Index);
                    h_pull = HistFac.makeTH1D(title.str(),
                                              "Pull", "#",
                                              pull_bins,
                                              "pull_"+varname);
                    pulls[varname] = h_pull;
                }
                else {
                    h_pull = it_pull->second;
                }
                h_pull->Fill(var.Pull);
            }
        }

        branches.b_fitok = (fitres.Status == APLCON::Result_Status_t::Success);

        branches.kinfit_chi2 = fitres.ChiSquare / fitres.NDoF;
        branches.b_fitIterations = unsigned(fitres.NIterations);

        branches.b_fitted_p = analysis::utils::ParticleVars(*fitter.GetFittedProton());


        TParticleList rec_photons(3);
        TParticlePtr  rec_proton = nullptr;
        TParticleList true_particles(4);

        branches.b_ggIM_real    = std::numeric_limits<double>::quiet_NaN();
        branches.b_ggIM_comb[0] = std::numeric_limits<double>::quiet_NaN();
        branches.b_ggIM_comb[1] = std::numeric_limits<double>::quiet_NaN();

        if(particletree && (branches.b_SigBgFlag==0 || branches.b_SigBgFlag ==1)) {
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

            branches.b_true_p = analysis::utils::ParticleVars(*true_particles[3]);

            const auto matched  = utils::match1to1(true_particles, data.Particles.GetAll(), [] (const TParticlePtr& p1, const TParticlePtr& p2) {
                return p1->Angle(p2->Vect());
            }, {0.0, degree_to_radian(15.0)});

            if(matched.size() == true_particles.size()) {

                rec_photons[0] = utils::FindMatched(matched, true_particles[0]);
                rec_photons[1] = utils::FindMatched(matched, true_particles[1]);
                rec_photons[2] = utils::FindMatched(matched, true_particles[2]);
                rec_proton     = utils::FindMatched(matched, true_particles[3]);

                FASSERT(rec_photons[0]!=nullptr);
                FASSERT(rec_photons[1]!=nullptr);
                FASSERT(rec_photons[2]!=nullptr);
                FASSERT(rec_proton    !=nullptr);


                if(rec_proton->Type() == ParticleTypeDatabase::Proton) {
                    branches.b_found_proton = 1;
                } else {
                    branches.b_found_proton = 0;
                }

                branches.b_ggIM_real    = IM(rec_photons[1], rec_photons[2]);
                branches.b_ggIM_comb[0] = IM(rec_photons[0], rec_photons[1]);
                branches.b_ggIM_comb[1] = IM(rec_photons[0], rec_photons[2]);

            }

        }


        tree->Fill();

    }

}

int OmegaEtaG2::identify(const TEvent &event)
{

    const auto particletree = event.MCTrue->ParticleTree;

    if(!particletree)
        return 0;

    const auto reaction_name = utils::ParticleTools::GetDecayString(particletree);

    for(const auto& c : reaction_channels) {

        if(particletree->IsEqual(c.second, utils::ParticleTools::MatchByParticleName)) {
            found_channels->Fill(reaction_name.c_str(), 1.0);
            return c.first;
        }
    }

    missed_channels->Fill(reaction_name.c_str(), 1.0);

    return 100000;
}

std::map<int, std::shared_ptr<OmegaEtaG2::decaytree_t> > OmegaEtaG2::makeChannels()
{
    std::map<int, std::shared_ptr<OmegaEtaG2::decaytree_t> > m;

    m[0] = ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Omega_gEta_3g);  //sig
    m[1] = ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Omega_gPi0_3g);  //ref

    m[2] = ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Direct1Pi0_2g);
    m[3] = ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Direct2Pi0_4g);
    m[4] = ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Direct3Pi0_6g);
    m[5] = ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::DirectPi0Eta_4g);
    m[6] = ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Omega_Pi0PiPPiM_2g);

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
    reaction_channels(makeChannels()),
    tree(HistFac.makeTTree("tree")),
    fitter("OmegaEtaG2", 3)
{
    const auto setup = ant::ExpConfig::Setup::GetLastFound();

    if(!setup) {
        throw std::runtime_error("No Setup found");
    }

    if(Options->Get<string>("Proton") == "MCTrue") {
        data_proton = false;
        LOG(INFO) << "Using proton from MCTrue";
    };

    if(Options->Get<string>("Tagger") == "MCTrue") {
        data_tagger = false;
        LOG(INFO) << "Using Tagger from MCTrue";
    };

    just_preselect = Options->Get<bool>("Preselect", false);

    if(just_preselect) {
        LOG(INFO) << "Running in preselect mode";
    }

    promptrandom.AddPromptRange({-5,5});
    promptrandom.AddRandomRange({-20, -10});
    promptrandom.AddRandomRange({ 10,  20});

    cut_ESum = Options->Get<double>("ESum", cut_ESum);

    branches.SetupBranches(tree);

    // set up kin fitter
    fitter.SetupBranches(tree, "EPB");

    fitter.LoadSigmaData(setup->GetPhysicsFilesDirectory()+"/FitterSigmas.root");

    signal_tree    = ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Omega_gEta_3g);
    reference_tree = ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Omega_gPi0_3g);

    steps = HistFac.makeTH1D("Steps","Step","Events passed",BinSettings(14),"steps");

    h_TotalEvents = HistFac.makeTH1D("TotalEvents","","",BinSettings(3),"TotalEvents");

    missed_channels = HistFac.makeTH1D("Unlisted Background Channels","","",BinSettings(20),"mBgC");
    found_channels  = HistFac.makeTH1D("Known Channels","","",BinSettings(reaction_channels.size()),"known");

    int i=1;
    for(const auto& c : reaction_channels) {
        found_channels->GetXaxis()->SetBinLabel(i++,utils::ParticleTools::GetDecayString(c.second).c_str());
    }

}

OmegaEtaG2::~OmegaEtaG2()
{

}

OmegaEtaG2::branches_t::branches_t() {}

void OmegaEtaG2::branches_t::SetupBranches(TTree* tree)
{
    b_g1.SetBranches(tree, "g1");
    b_g2.SetBranches(tree, "g2");
    b_g3.SetBranches(tree, "g3");
    b_p.SetBranches(tree, "p");
    b_ggg.SetBranches(tree, "ggg");
    b_mmvector.SetBranches(tree,"mmvect");
    b_fitted_p.SetBranches(tree, "fitted_p");
    b_true_p.SetBranches(tree, "ture_p");

    tree->Branch("p_Time",      &b_pTime);
    tree->Branch("p_PSA_R",     &b_p_PSA_R);
    tree->Branch("p_PSA_Angle", &b_p_PSA_Angle);
    tree->Branch("p_detector",  &b_p_detector);

    tree->Branch("ggg_Time", &b_gggTime);
    tree->Branch("ggIM",     b_ggIM, "ggIM[3]/D");

    tree->Branch("CoplAngle", &b_copl_angle);
    tree->Branch("Angle_p_mm",&b_p_mm_angle);

    tree->Branch("found_proton", &b_found_proton);

    tree->Branch("TagCh",   &b_TagCh);
    tree->Branch("TagE",    &b_TagE);
    tree->Branch("TagTime", &b_TagTime);
    tree->Branch("TagW",    &b_TagW);

    tree->Branch("SigBgFlag",       &b_SigBgFlag);

    tree->Branch("BachelorE[3]",     b_BachelorE, "BachelorE[3]/D");

    tree->Branch("CBAvgTime",       &b_CBAvgTime);

    tree->Branch("ggIM_real",       &b_ggIM_real);
    tree->Branch("ggIM_comb[2]",     b_ggIM_comb, "ggIM_comb[2]/D");

    tree->Branch("chi2dof",         &kinfit_chi2);
    tree->Branch("fitstatus",       &b_fitok);
    tree->Branch("fitIterations",   &b_fitIterations);


}


AUTO_REGISTER_PHYSICS(OmegaEtaG)
AUTO_REGISTER_PHYSICS(OmegaMCTruePlots)
AUTO_REGISTER_PHYSICS(OmegaMCTree)
AUTO_REGISTER_PHYSICS(OmegaEtaG2)
