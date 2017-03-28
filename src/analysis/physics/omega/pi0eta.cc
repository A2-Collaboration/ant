#include "pi0eta.h"
#include "TH1D.h"
#include "utils/Combinatorics.h"
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

#include "utils/ParticleTools.h"
#include "utils/Matcher.h"

#include "APLCON.hpp"
#include "expconfig/ExpConfig.h"
#include "base/WrapTFile.h"
#include "TCanvas.h"
#include <cassert>

#include "utils/Matcher.h"

#include "utils/uncertainties/Optimized.h"

#include "root-addons/analysis_codes/hstack.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace ant::std_ext;


void Pi0Eta::Analyse(const TEventData &data, const TEvent& event, manager_t& manager)
{

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

    if(!triggersimu.HasTriggered())
        return;

    steps->Fill("1 Triggered", 1);

    const auto n_cands = geoAccepted(data.Candidates);

    if(n_cands != 5) {
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

    if(iphotons.size() != 4)
        return;

    if(iprotons.size() != 1)
        return;

    const TParticleList protons = FilterProtons(getGeoAccepted(iprotons));

    if(protons.size() != 1)
        return;

    const TParticleList photons = FilterPhotons(getGeoAccepted(iphotons));

    if(photons.size() != 4)
        return;

    steps->Fill("3 nPhotons nProtons", 1);

    const auto& proton = protons.at(0);

    t.photons().at(0) = *photons.at(0);
    t.photons().at(1) = *photons.at(1);
    t.photons().at(2) = *photons.at(2);

    t.p      = *proton;
    t.p_Time = getTime(proton);

    t.p_detector   = 0;

    if(proton->Candidate) {

        if(proton->Candidate->Detector & Detector_t::Type_t::TAPS) {
            t.p_detector = 2;
        } else if(proton->Candidate->Detector & Detector_t::Type_t::CB) {
            t.p_detector = 1;
        }

    }


    const TParticle gggg(ParticleTypeDatabase::Omega, LVSum(photons.begin(), photons.end()));
    t.gggg = gggg;

    t.copl_angle = fabs(vec2::Phi_mpi_pi(proton->Phi() - gggg.Phi() - M_PI));

    if(t.copl_angle > cut_Copl)
        return;

    steps->Fill("4 Coplanarity", 1);

    t.CBAvgTime = triggersimu.GetRefTiming();
    if(!isfinite(t.CBAvgTime))
        return;

    steps->Fill("5 valid CB Avg Time", 1);

    if(data.TaggerHits.size() > 0)
        steps->Fill("6 has TaggHits", 1);

    for(const TTaggerHit& TagH : data.TaggerHits) {

        promptrandom.SetTaggerTime(triggersimu.GetCorrectedTaggerTime(TagH));

        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        t.TaggW  = promptrandom.FillWeight();
        t.TaggE  = TagH.PhotonEnergy;
        t.TaggCh = TagH.Channel;
        t.TaggT  = TagH.Time;

        const LorentzVec beam_target = TagH.GetPhotonBeam() + LorentzVec({0, 0, 0}, ParticleTypeDatabase::Proton.Mass()); // make global
        const TParticle missing(ParticleTypeDatabase::Proton, beam_target - gggg);

        t.mm = missing;

        t.p_mm_angle = radian_to_degree(missing.Angle(*proton));

        TParticleList kinfitted_photons;
        // KinFit
        {
            auto fitres = fitter.DoFit(TagH.PhotonEnergy, proton, photons);

            if(fitres.Status != APLCON::Result_Status_t::Success)
                continue;

            t.KinFitChi2 = fitres.ChiSquare / fitres.NDoF;
            t.KinFitProb = fitres.Probability;
            t.KinFitIterations = unsigned(fitres.NIterations);

            t.p_fitted = *fitter.GetFittedProton();

            kinfitted_photons = fitter.GetFittedPhotons();

            t.photons_fitted().at(0) = *kinfitted_photons.at(0);
            t.photons_fitted().at(1) = *kinfitted_photons.at(1);
            t.photons_fitted().at(2) = *kinfitted_photons.at(2);
            t.photons_fitted().at(3) = *kinfitted_photons.at(3);


        }



        const TParticle gggg_fitted(ParticleTypeDatabase::Omega, LVSumL(t.photons_fitted().begin(), t.photons_fitted().end()));
        t.gggg_fitted = gggg_fitted;

        //===== Hypothesis testing with kinematic fitter ======

        {


            // Kin fit: test pi0 hypothesis

            treefitter.treefitter.PrepareFits(TagH.PhotonEnergy, proton, photons);

            treefitter.HypTestCombis(photons, kinfitted_photons,
                                     t.TreeFitChi2,
                                     t.TreeFitProb,
                                     t.TreeFitpi0,
                                     t.TreeFiteta,
                                     t.TreeFitgggg);

            if(opt_save_after_treefit && (t.TreeFitChi2 < opt_treefit_chi2cut)) {
                manager.SaveEvent();
            }


        }


        tree->Fill();

    }

}


Pi0Eta::ReactionChannelList_t Pi0Eta::makeChannels()
{
    ReactionChannelList_t m;

    m.channels[0] = ReactionChannel_t("Data");
    m.channels[1] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0Eta_4g), kRed};  //sig
    m.channels[2] = ReactionChannel_t("Sum MC");
    m.channels[3] = ReactionChannel_t("MC BackG");

    m.channels[10] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0_2g),      "#pi^{0}", kYellow};
    m.channels[11] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::TwoPi0_4g),   "#pi^{0} #pi^{0}", kOrange};
    m.channels[12] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::ThreePi0_6g), "#pi^{0} #pi^{0} #pi^{0}",kGreen-9};
    m.channels[13] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0Eta_4g),   "#pi^{0} #eta",kBlue};
    m.channels[14] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Omega_Pi0PiPPiM_2g),"#omega #rightarrow #pi^{0} #pi^{+} #pi^{-}",kMagenta};
    m.channels[m.other_index] = ReactionChannel_t(nullptr, "Others", kCyan);

    return m;
}

bool Pi0Eta::AcceptedPhoton(const TParticlePtr& photon)
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

bool Pi0Eta::AcceptedProton(const TParticlePtr& proton)
{
    if(proton_theta.Contains(proton->Theta())){
        return true;
    }

    return false;
}

TParticleList Pi0Eta::FilterPhotons(const TParticleList& list)
{

    TParticleList olist;

    for(const auto& p : list) {
        if(AcceptedPhoton(p)) {
            olist.emplace_back(p);
        }
    }
    return olist;
}

TParticleList Pi0Eta::FilterProtons(const TParticleList& list)
{

    TParticleList olist;

    for(const auto& p : list) {
        if(AcceptedProton(p)) {
            olist.emplace_back(p);
        }
    }
    return olist;
}


Pi0Eta::Pi0Eta(const std::string& name, OptionsPtr opts):
    OmegaBase(name, opts),
    tree(HistFac.makeTTree("tree")),

    cut_ESum(opts->Get<double>("CBESum", 550.0)),
    cut_Copl(degree_to_radian(opts->Get<double>("CoplAngle", 15.0))),
    photon_E_cb(opts->Get<decltype(photon_E_cb)>("PhotonECB", {50.0, 1600.0})),
    photon_E_taps(opts->Get<decltype(photon_E_taps)>("PhotonETAPS", {200.0, 1600.0})),
    proton_theta(degree_to_radian(opts->Get<decltype(proton_theta)>("ProtonThetaRange", {2.0, 45.0}))),
    model(make_shared<utils::UncertaintyModels::Optimized_Oli1>()),
    fitter(model),
    treefitter(
        ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0Eta_4g),
        model
        ),
    opt_save_after_treefit(opts->Get("SaveAfterTreefit", false)),
    opt_treefit_chi2cut(opts->Get<double>("TreeFit_Chi2Cut", 10.0))
{

    promptrandom.AddPromptRange({-5,5});
    promptrandom.AddRandomRange({-20, -10});
    promptrandom.AddRandomRange({ 10,  20});

    t.CreateBranches(tree);

    missed_channels = HistFac.makeTH1D("Unlisted Channels", "", "Total Events seen", BinSettings(20), "unlistedChannels");
    found_channels  = HistFac.makeTH1D("Listed Channels",   "", "Total Events seen", BinSettings(20), "listedChannels");

    for(const auto& c : reaction_channels.channels) {

        stephists[c.first] = HistFac.makeTH1D("Steps: " + c.second.name, "", "", BinSettings(14), "steps_" + to_string(c.first));
        stephists[c.first]->SetLineColor(c.second.color);

        if(c.first<20)
            found_channels->GetXaxis()->SetBinLabel(c.first+1,c.second.name.c_str());
    }

    //fitter.SetupBranches(t.Tree);

}

Pi0Eta::~Pi0Eta()
{
}

void Pi0Eta::Finish()
{
}



Pi0Eta::Pi0EtaTree_t::Pi0EtaTree_t()
{}



Pi0Eta::ReactionChannel_t::ReactionChannel_t(const std::shared_ptr<OmegaEtaG2::decaytree_t> &t, const string &n, const int c):
    name(n), tree(t), color(c)
{
}

Pi0Eta::ReactionChannel_t::ReactionChannel_t(const std::shared_ptr<OmegaEtaG2::decaytree_t> &t, const int c):
    name(utils::ParticleTools::GetDecayString(t)),
    tree(t),
    color(c)
{}

Pi0Eta::ReactionChannel_t::ReactionChannel_t(const string &n):
    name(n)
{}

Pi0Eta::ReactionChannel_t::~ReactionChannel_t()
{}



unsigned Pi0Eta::ReactionChannelList_t::identify(const ant::TParticleTree_t& tree) const
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

const Pi0Eta::ReactionChannelList_t Pi0Eta::reaction_channels = Pi0Eta::makeChannels();

const unsigned Pi0Eta::ReactionChannelList_t::other_index = 1000;



Pi0Eta::MyTreeFitter_t::MyTreeFitter_t(const ParticleTypeTree& ttree, utils::UncertaintyModelPtr model):
    treefitter(
        ttree,
        model, false,
        [] (const ParticleTypeTree& t) { return utils::TreeFitter::nodesetup_t((t->Get() == ParticleTypeDatabase::Omega)); }
        )
{

    fitted_pi0     = treefitter.GetTreeNode( ParticleTypeDatabase::Pi0);
    fitted_pi0_g1  = fitted_pi0->Daughters().front();
    fitted_pi0_g2  = fitted_pi0->Daughters().back();

    fitted_eta     = treefitter.GetTreeNode(ParticleTypeDatabase::Eta);
    fitted_eta_g1  = fitted_eta->Daughters().front();
    fitted_eta_g2  = fitted_eta->Daughters().back();

    if(!fitted_pi0 || !fitted_eta)
        throw std::runtime_error("Error initializing Pi0Eta::MyTreeFitter_t");
}

void Pi0Eta::MyTreeFitter_t::HypTestCombis(const TParticleList& unfitted, const TParticleList& kinfitted,
                           double& chi2,
                           double& prob,
                           TLorentzVector& pi0,
                           TLorentzVector& eta,
                           TLorentzVector& gggg)
{

    auto getIndex = [] (const TParticlePtr& p, const TParticleList& plist) {
        for(size_t i=0; i<plist.size(); ++i ) {
            if(plist.at(i) == p)
                return i;
        }

        throw std::runtime_error("Partile not found in list");
    };

    APLCON::Result_t treefitres;

    chi2 = inf;
    prob = -inf;

    while(treefitter.NextFit(treefitres)) {

        const auto _chi2 = treefitres.Status == APLCON::Result_Status_t::Success ? treefitres.ChiSquare : NaN;

        if(isfinite(_chi2) && (_chi2 < chi2)) {
            chi2 = _chi2;
            prob = treefitres.Status == APLCON::Result_Status_t::Success ? treefitres.Probability : NaN;

            {
                const TParticlePtr& g1 = fitted_pi0_g1->Get().Leaf->Particle;
                const TParticlePtr& g2 = fitted_pi0_g2->Get().Leaf->Particle;
                pi0   = *kinfitted.at(getIndex(g1,unfitted)) +  *kinfitted.at(getIndex(g2,unfitted));
            }

            {
                const TParticlePtr& g1 = fitted_eta_g1->Get().Leaf->Particle;
                const TParticlePtr& g2 = fitted_eta_g2->Get().Leaf->Particle;
                eta   = *kinfitted.at(getIndex(g1,unfitted)) +  *kinfitted.at(getIndex(g2,unfitted));
            }
            gggg = pi0 + eta;
        }
    }
}

Pi0Eta::MyTreeFitter_t::~MyTreeFitter_t()
{

}

AUTO_REGISTER_PHYSICS(Pi0Eta)
