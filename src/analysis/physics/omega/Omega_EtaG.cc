#include "Omega_EtaG.h"
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

#include "analysis/utils/uncertainties/Interpolated.h"
#include "analysis/utils/uncertainties/FitterSergey.h"
#include "analysis/utils/uncertainties/Optimized.h"
#include "analysis/utils/uncertainties/MeasuredProton.h"
#include "analysis/plot/CutTree.h"
#include "base/vec/vec2.h"
#include "TStyle.h"
#include "TCutG.h"

#include "root-addons/analysis_codes/hstack.h"

#include "analysis/physics/Plotter.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace ant::std_ext;

void OmegaBase::ProcessEvent(const TEvent& event, manager_t& manager)
{
    triggersimu.ProcessEvent(event);
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

string to_string(const OmegaBase::DataMode &m)
{
    if(m == OmegaBase::DataMode::MCTrue) {
        return "MCTrue";
    } else {
        return "Reconstructed";
    }
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

#define FASSERT(x) if(!(x)) LOG(ERROR) << "ERROR";


double IM(const TParticlePtr& p1, const TParticlePtr& p2) {
    return (*p1+*p2).M();
}

double getTime(const TParticlePtr& p) {
    return p->Candidate != nullptr ? p->Candidate->Time : std_ext::NaN;
}



/**
 * @brief copy TCandidateList to new TCandidatePtrList
 * @param c
 * @return
 */
TCandidatePtrList copy(const TCandidateList& c) {

    TCandidatePtrList res(c.size());
    copy(c.get_iter().begin(), c.get_iter().end(), res.begin());
    return res;
}

void SortByEnergy(TCandidatePtrList& c) {
    sort(c.begin(), c.end(), [] (const TCandidatePtr& a, const TCandidatePtr& b) { return a->CaloEnergy > b->CaloEnergy; });
}


bool OmegaEtaG2::ProtonCheck(const TCandidatePtr& c) const
{
    return proton_theta.Contains(c->Theta);
}

bool OmegaEtaG2::PhotonCheck(const TCandidatePtr& c) const {

    if(c->Detector & Detector_t::Type_t::CB)
        return photon_E_cb.Contains(c->CaloEnergy);

    if(c->Detector & Detector_t::Type_t::TAPS)
        return photon_E_taps.Contains(c->CaloEnergy);

    return false;
}


static TVector2 getPSAVector(const TParticlePtr& p) {
    if(p->Candidate) {
        const auto cluster = p->Candidate->FindCaloCluster();
        if(cluster) {
            return {cluster->Energy, cluster->ShortEnergy};
        }
    }

    throw std::runtime_error("Incomplete Particle without candiate or CaloCluster");

}

static int getDetectorAsInt(const Detector_t::Any_t& d) {

    if(d & Detector_t::Type_t::TAPS) {
        return 2;
    } else if(d & Detector_t::Type_t::CB) {
        return 1;
    }

    return 0;
}

/**
 * @brief StrictPhotonVeto
 * @param photon
 * @param proton
 *
 * in TAPS -> accept
 * in CB:
 *   uncharged -> accept
 *   charged:
 *      same PID element as proton -> accept
 *      different element -> reject
 * @return
 */
bool OmegaEtaG2::StrictPhotonVeto(const TCandidate& photon, const TCandidate& proton) const {

    if(photon.Detector & Detector_t::Type_t::CB) {

        if(photon.VetoEnergy < 0.25)
            return true;

        if(fabs(vec2::Phi_mpi_pi(proton.Phi - photon.Phi - M_PI)) < degree_to_radian(15.0)) {
            return true;
        }

        return false;
    }

    return true;
}

struct MaxTracker {
    double v;
    MaxTracker(const double& start = -std_ext::inf) : v(start) {}
    bool Track(const double& value) { if(value > v) {v = value; return true;} else return false; }
    double operator()() const { return v; }
};

void OmegaEtaG2::Analyse(const TEventData &data, const TEvent& event, manager_t&) {

    dCounters.EventStart();
    const bool is_MC = event.Reconstructed().ID.isSet(TID::Flags_t::MC);
    const auto& mctree = event.MCTrue().ParticleTree;

    t.Channel = reaction_channels.identify(mctree);
    t.ChannelString = utils::ParticleTools::GetDecayString(mctree);

    if(t.Channel == ReactionChannelList_t::other_index) {
        if(event.MCTrue().ParticleTree!=nullptr) {
            missed_channels->Fill(t.ChannelString().c_str(), 1.0);
        }
    } else {
        found_channels->Fill(t.Channel);
    }

    TH1* steps = stephists.at(t.Channel);

    steps->Fill("Events seen", 1);


    if(data.Candidates.size() < nCandsMin)
        return;

    steps->Fill("NCans OK", 1);

    t.nCandsInput = unsigned(data.Candidates.size());

    TCandidatePtrList cands;
    for(const auto& c : data.Candidates.get_iter()) {
        if(cand_theta_range.Contains(c->Theta))
            cands.emplace_back(c);
    }

    if(cands.size() < nCandsMin || cands.size() > nCandsMax)
        return;

    steps->Fill("NCans after cut OK", 1);

    SortByEnergy(cands);

    t.CandsUsedE   = 0.0;
    t.CandsunUsedE = 0.0;
    t.nTouchesHole = 0;

    for(size_t i=0; i<cands.size(); ++i) {

        const auto& c = cands.at(i);

        if(i<nCandsMin)
            t.CandsUsedE   += c->CaloEnergy;
        else
            t.CandsunUsedE += c->CaloEnergy;

        if(c->FindCaloCluster()->HasFlag(TCluster::Flags_t::TouchesHoleCentral))
            t.nTouchesHole += 1;
    }

    cands.resize(nCandsMin);

    t.CBESum = triggersimu.GetCBEnergySum();

    if(t.CBESum < cut_ESum)
        return;
    steps->Fill("CB ESum OK", 1);

    t.CBAvgTime = triggersimu.GetRefTiming();
    if(!isfinite(t.CBAvgTime))
        return;

    steps->Fill("CB Avg Time OK", 1);



    if(data.TaggerHits.size() > 0)
        steps->Fill("6 has TaggHits", 1);


    {
        const auto mctrue_particles = utils::ParticleTypeList::Make(mctree);
        const auto& mctrue_protons = mctrue_particles.Get(ParticleTypeDatabase::Proton);
        if(mctrue_protons.size() == 1) {
            t.p_true = *(mctrue_protons.front());
        } else {
            t.p_true = {};
        }
    }

    tagChMult.Fill(data.TaggerHits);

    auto dcTaggerHitsAccepted = dTaggerHitsAccepted.getHandle();

    // set fitter uncertainty models
    {
        const auto& model = is_MC ? model_MC : model_Data;
        fitter.SetUncertaintyModel(model);
        fitter_pi0.treefitter.SetUncertaintyModel(model);
        fitter_eta.treefitter.SetUncertaintyModel(model);
    }

    for(const TTaggerHit& TagH : data.TaggerHits) {

        dCounters.TaggerLoopBegin();

        promptrandom.SetTaggerTime(triggersimu.GetCorrectedTaggerTime(TagH));

        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        dCounters.TaggerHitAccepted();

        dcTaggerHitsAccepted.Count();


        t.TaggW  = promptrandom.FillWeight();
        t.TaggE  = TagH.PhotonEnergy;
        t.TaggCh = TagH.Channel;
        t.TaggT  = TagH.Time;


        MaxTracker kinfit_best_prob(opt_kinfit_probcut);

        TParticleList selected_photons;
        TParticleList fitted_photons;
        TParticlePtr  selected_proton;
        TParticlePtr  fitted_proton;
        bool fit_ok = false;

        for(auto it_proton = cands.begin(); it_proton != cands.end(); ++it_proton) {

            dCounters.PIDLoopBegin();

            //proton acceptance checks
            if(!ProtonCheck(*it_proton))
                continue; //proton loop

            const TParticlePtr proton = make_shared<TParticle>(ParticleTypeDatabase::Proton, *it_proton);


            TParticleList photons;
            photons.reserve(nphotons);
            for(auto it_photon = cands.begin(); it_photon!=cands.end(); ++it_photon) {

                if(it_photon == it_proton)
                    continue;

                if(PhotonCheck(*it_proton)) {

                    if(!opt_strict_Vetos || StrictPhotonVeto(**it_photon, **it_proton)) {
                        photons.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::Photon, *it_photon));
                    }

                }
            }

            if(photons.size() != nphotons)
                continue; //proton loop

            dCounters.ParticlesFound();


            const TParticle ggg(ParticleTypeDatabase::Omega, LVSum(photons.begin(), photons.end()));

            const LorentzVec beam_target = TagH.GetPhotonBeam() + target;
            const TParticle missing_vector(ParticleTypeDatabase::Proton, beam_target - ggg);

            const auto p_mm_angle = proton->Angle(missing_vector);

            if(!cut_missing_mass.Contains(missing_vector.M()))
                continue; //proton loop

            dCounters.MissingMassOK();

            dCounters.KinfitLoopBegin();
            const auto fitres = fitter.DoFit(TagH.PhotonEnergy, proton, photons);

            if(fitres.Status != APLCON::Result_Status_t::Success)
                continue; //proton loop
            steps->Fill("Fit success", 1.0);
            hKinfitProb_all->Fill(fitres.Probability);

            if(kinfit_best_prob.Track(fitres.Probability)) {

                dCounters.KinfitHighscore();

                fit_ok = true;

                // update Tree branches with new best values

                // Kin Fit
                t.KinFitChi2       = fitres.ChiSquare / fitres.NDoF;
                t.KinFitProb       = fitres.Probability;
                t.KinFitIterations = fitres.NIterations;

                // proton
                selected_proton = proton;
                fitted_proton   = fitter.GetFittedProton();

                // photons
                fitted_photons   = fitter.GetFittedPhotons();
                selected_photons = photons;

                t.p          = *selected_proton;
                t.p_fitted   = *fitted_proton;
                t.p_Time     = selected_proton->Candidate->Time;
                t.p_PSA      = getPSAVector(selected_proton);
                t.p_vetoE    = selected_proton->Candidate->VetoEnergy;
                t.p_detector = getDetectorAsInt(selected_proton->Candidate->Detector);

                const auto fitparticles = fitter.GetFitParticles();
                assert(fitparticles.size() == nphotons +1);

                t.p_theta_pull  = fitparticles.at(0).GetPulls().at(1);
                t.p_phi_pull    = fitparticles.at(0).GetPulls().at(2);

                for(size_t i=0; i<nphotons; ++i) {
                    t.photons().at(i)            = *(selected_photons.at(i));
                    t.photons_fitted().at(i)     = *(fitted_photons.at(i));
                    t.photons_Time().at(i)       = selected_photons.at(i)->Candidate->Time;
                    t.photons_vetoE().at(i)      = selected_photons.at(i)->Candidate->VetoEnergy;
                    t.photons_PSA().at(i)        = getPSAVector(selected_photons.at(i));
                    t.photons_detector().at(i)   = getDetectorAsInt(selected_photons.at(i)->Candidate->Detector);
                    t.photon_E_pulls().at(i)     = fitparticles.at(i+1).GetPulls().at(0);
                    t.photon_theta_pulls().at(i) = fitparticles.at(i+1).GetPulls().at(1);
                    t.photon_phi_pulls().at(i)   = fitparticles.at(i+1).GetPulls().at(2);
                }

                // other
                t.beam_E_pull = fitter.GetBeamEPull();
                t.beam_E_fitted = fitter.GetFittedBeamE();

                t.zVertex = fitter.GetFittedZVertex();

                t.ggg() = ggg;

                const TParticle ggg_fitted(ParticleTypeDatabase::Omega, LVSum(fitted_photons.begin(), fitted_photons.end()));
                t.ggg_fitted = ggg_fitted;

                t.mm()  = missing_vector;
                t.p_mm_angle = p_mm_angle;
                t.copl_angle = fabs( vec2::Phi_mpi_pi(proton->Phi() - ggg.Phi() - M_PI));

                const auto gggBoost = -ggg.BoostVector();
                const auto gggBoost_fitted = -ggg_fitted.BoostVector();


                const auto bt = LorentzVec(TagH.GetPhotonBeam().p,TagH.GetPhotonBeam().E +ParticleTypeDatabase::Proton.Mass());

                const auto cm_boost = bt.BoostVector();

                LorentzVec ggg_cm;
                ggg_cm.Boost(-cm_boost);

                t.ggg_cm = ggg_cm;

                for(const auto& comb : combs) {
                    const auto& combindex = comb[2];

                    {
                        const auto& g1 = photons.at(comb[0]);
                        const auto& g2 = photons.at(comb[1]);
                        const auto& g3 = photons.at(comb[2]);

                        const auto gg = *g1 + *g2;
                        t.ggIM().at(combindex) = gg;
                        const auto g3_boosted = Boost(*g3, gggBoost);
                        t.BachelorE().at(combindex) = g3_boosted.E;
                    }

                    {
                        const auto& g1 = t.photons_fitted().at(comb[0]);
                        const auto& g2 = t.photons_fitted().at(comb[1]);
                        const auto& g3 = t.photons_fitted().at(comb[2]);

                        const auto gg = g1 + g2;

                        t.ggIM_fitted().at(combindex) = gg;

                        const auto g3_boosted = Boost(g3, gggBoost_fitted);

                        t.BachelorE_fitted().at(combindex) = g3_boosted.E;
                    }

                }
            }
        } // Proton Loop

        if(!fit_ok)
            continue;

        dCounters.AfterKinfitLoop();

        steps->Fill("Kinfit OK", 1.0);
        if(!opt_skip_treefit) {

        //===== Hypothesis testing with kinematic fitter ======

        {

            // Kin fit: test pi0 hypothesis
            dCounters.TreeFit();

            fitter_pi0.treefitter.PrepareFits(TagH.PhotonEnergy, selected_proton, selected_photons);

            fitter_pi0.HypTestCombis(selected_photons,
                                     t.pi0chi2,
                                     t.pi0prob,
                                     t.pi0_im,
                                     t.pi0_omega_im,
                                     t.iBestPi0);



            // Kin fit: test eta hypothesis
            dCounters.TreeFit();

            fitter_eta.treefitter.PrepareFits(TagH.PhotonEnergy, selected_proton, selected_photons);

            fitter_eta.HypTestCombis(selected_photons,
                                     t.etachi2,
                                     t.etaprob,
                                     t.eta_im,
                                     t.eta_omega_im,
                                     t.iBestEta);


            // find most probable hypothesis

            t.bestHyp = 0;

            // both fits failed
            if(t.iBestEta == -1 && t.iBestPi0 == -1)
                continue;

            // Pi0 fits failed, but at least one eta fit worked -> eta wins
            if(t.iBestPi0 == -1) {
                t.bestHyp = 2; // Eta
            }

            // Eta fits failed, but at least one pi0 fit worked -> pi0 wins
            else if(t.iBestEta == -1) {
                t.bestHyp = 1; // Pi0
            }

            // both hypotheses have valid fit results: select lower chi2
            else {
                if(t.pi0chi2().at(t.iBestPi0) < t.etachi2().at(t.iBestEta)) {
                    t.bestHyp = 1;  // PI0
                } else {
                    t.bestHyp = 2;  // ETA
                }
            }

        }
        }

        dCounters.TaggerLoopEnd();
        tree->Fill();

    } // Tagger Loop

    dCounters.EventEnd();

}

size_t OmegaEtaG2::CombIndex(const TParticleList& orig, const MyTreeFitter_t& f)
{
    for(size_t i=0; i<orig.size(); ++i) {
        if(orig[i] == f.fitted_g_Omega->Get().Leaf->Particle) {
            return i;
        }
    }

    throw std::runtime_error("CombIndex: Photon not found");

}

OmegaEtaG2::ReactionChannelList_t OmegaEtaG2::makeChannels()
{
    ReactionChannelList_t m;

    m.channels[0] = ReactionChannel_t("Data");
    m.channels[1] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Omega_gEta_3g), kRed};  //sig
    m.channels[2] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Omega_gPi0_3g), kGreen};  //ref
    m.channels[3] = ReactionChannel_t("Sum MC");
    m.channels[4] = ReactionChannel_t("MC BackG");


    m.channels[10] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0_2g),            "#pi^{0}",                                        kYellow};
    m.channels[11] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::TwoPi0_4g),         "#pi^{0} #pi^{0}",                                kOrange};
    m.channels[12] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::ThreePi0_6g),       "#pi^{0} #pi^{0} #pi^{0}",                        kGreen-9};
    m.channels[13] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0Eta_4g),         "#pi^{0} #eta",                                   kBlue};
    m.channels[14] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Omega_Pi0PiPPiM_2g),"#omega #rightarrow #pi^{0} #pi^{+} #pi^{-}",     kMagenta};
    m.channels[15] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Rho_PiPi),          "#rho #rightarrow #pi^{+} #pi^{-}",               kTeal};
    m.channels[16] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0PiPi_2gPiPi),    "#pi^{0} (#rightarrow 2 #gamma) #pi^{+} #pi^{-}", kPink};
    m.channels[17] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::TwoPi0PiPi_4gPiPi), "#pi^{0} #pi^{0} #pi^{+} #pi^{-}",                kRed};
    m.channels[18] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Eta_2g),            "#eta #rightarrow 2#gamma",                       kGreen};
    m.channels[m.other_index] = ReactionChannel_t(nullptr, "Others", kCyan);

    return m;
}

OmegaEtaG2::OmegaEtaG2(const std::string& name, OptionsPtr opts):
    OmegaBase(opts->Get<string>("Name", name), opts),
    tree(HistFac.makeTTree("tree")),
    cut_ESum(                          opts->Get<double>(                    "CBESum",                    600.0)),
    cut_Angle_PMM(degree_to_radian(    opts->Get<double>(                    "PMM_angle",             20.0))),
    photon_E_cb(                       opts->Get<decltype(photon_E_cb)>  (   "PhotonECB",        { 0.0,  1600.0})),
    photon_E_taps(                     opts->Get<decltype(photon_E_taps)>(   "PhotonETAPS",      { 50.0, 1600.0})),
    proton_theta(degree_to_radian(     opts->Get<decltype(proton_theta)> (   "ProtonThetaRange", { 7.0,    45.0}))),
    cand_theta_range(degree_to_radian( opts->Get<decltype(proton_theta)> (   "CandThetaRange",   { 7.0,   160.0}))),
    cut_missing_mass(                  opts->Get<decltype(cut_missing_mass)>("MissingMassWindow", interval<double>::CenterWidth(ParticleTypeDatabase::Proton.Mass(), 450.0))),
    cut_gggim(                         opts->Get<decltype(cut_gggim)>(       "GGGim",             interval<double>(500.0, std_ext::inf))),
    opt_kinfit_probcut(                opts->Get<double>(                    "KinFit_ProbCUt",       0.005)),
    opt_FitZVertex(                    opts->Get<bool>(                      "KinFit_FitVertex",     true)),
    opt_strict_Vetos(                  opts->Get<bool>(                      "Strict_Vetos",         true)),
    opt_z_sigma(                       opts->Get<double>(                    "ZSigma",               3.0)),
    opt_skip_treefit(		opts->Get<bool>( "SkipTreefit", false)),


    promptrandom(ExpConfig::Setup::Get()),
    model_Data(utils::UncertaintyModels::Interpolated::makeAndLoad(
                   utils::UncertaintyModels::Interpolated::Type_t::Data
                   )),
    model_MC(utils::UncertaintyModels::Interpolated::makeAndLoad(
                   utils::UncertaintyModels::Interpolated::Type_t::MC
                   )),
    fitter(model_Data, opt_FitZVertex),
    fitter_pi0(
        ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Omega_gPi0_3g),
        ParticleTypeDatabase::Pi0,
        model_Data,
        opt_FitZVertex,
        opts->Get<bool>(                      "KinFit_FixOmegaMass",     false)
        ),
    fitter_eta(
        ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Omega_gEta_3g),
        ParticleTypeDatabase::Eta,
        model_Data,
        opt_FitZVertex,
        opts->Get<bool>(                      "KinFit_FixOmegaMass",     false)
        ),
    tagChMult(HistFac),
    dTaggerHitsAccepted(HistFac.makeTH1D("Tagger Hits Accepted Per Event","","",BinSettings(10),"dTHAcceptedperEvent")),
    dCounters(HistFac)

{
    // initbialize fitter Z Vertex sigma
    {
        if(fitter.IsZVertexFitEnabled())
            fitter.SetZVertexSigma(opt_z_sigma);
        if(fitter_eta.treefitter.IsZVertexFitEnabled())
            fitter_eta.treefitter.SetZVertexSigma(opt_z_sigma);
        if(fitter_pi0.treefitter.IsZVertexFitEnabled())
            fitter_pi0.treefitter.SetZVertexSigma(opt_z_sigma);
    }

    t.CreateBranches(tree);

    missed_channels = HistFac.makeTH1D("Unlisted Channels","","Total Events seen",BinSettings(20),"unlistedChannels");
    found_channels  = HistFac.makeTH1D("Listed Channels",  "","Total Events seen",BinSettings(20),"listedChannels");

    for(const auto& c : reaction_channels.channels) {

        stephists[c.first] = HistFac.makeTH1D("Steps: " + c.second.name, "", "", BinSettings(14), "steps_" + to_string(c.first));
        stephists[c.first]->SetLineColor(short(c.second.color));

        if(c.first<20)
            found_channels->GetXaxis()->SetBinLabel(int(c.first+1),c.second.name.c_str());
    }

    hKinfitProb_all = HistFac.makeTH1D("Kinfit Prob", "Prob", "", BinSettings(1000,0,1),"kfprob");

    LOG(INFO) << "Initialized " << GetName() << ":";
    LOG(INFO) << " CBESum Cut           " << cut_ESum                   << " MeV";
    LOG(INFO) << " Coplanarity Cut      " << radian_to_degree(cut_Angle_PMM) << " deg";
    LOG(INFO) << " ggg IM cut           " << cut_gggim << " MeV";
    LOG(INFO) << " Min. Photon E (CB)   " << photon_E_cb                << " MeV";
    LOG(INFO) << " Min. Photon E (TAPS) " << photon_E_taps              << " MeV";
    LOG(INFO) << " Proton Theta Angle   " << radian_to_degree(proton_theta) << " deg";
    LOG(INFO) << " Missing Mass Window  " << cut_missing_mass           << " MeV";
    LOG(INFO) << " Min. Kinfit Prob     " << opt_kinfit_probcut;
    LOG(INFO) << " Strict Vetos         " << (opt_strict_Vetos ? "On" : "Off");
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

utils::TreeFitter::nodesetup_t::getter fgetter(const bool& fixmass) {

    if(fixmass)
        return {};

    return [] (const ParticleTypeTree& t) { return utils::TreeFitter::nodesetup_t((t->Get() == ParticleTypeDatabase::Omega)); };

}

OmegaEtaG2::MyTreeFitter_t::MyTreeFitter_t(const ParticleTypeTree& ttree, const ParticleTypeDatabase::Type& mesonT, utils::UncertaintyModelPtr model, const bool fix_z_vertex, const bool fix_omega_mass):
    treefitter(
        ttree,
        model, fix_z_vertex,
        fgetter(fix_omega_mass)
        )
{
    auto find_daughter = [] (utils::TreeFitter::tree_t& node, const ParticleTypeDatabase::Type& type) {
        for(const auto& d : node->Daughters()) {
            if(d->Get().TypeTree->Get() == type)
                return d;
        }
        return utils::TreeFitter::tree_t(nullptr);
    };

    fitted_Omega   = treefitter.GetTreeNode(ParticleTypeDatabase::Omega);
    fitted_g_Omega = find_daughter(fitted_Omega, ParticleTypeDatabase::Photon);
    fitted_X       = find_daughter(fitted_Omega, mesonT);
    fitted_g1_X    = fitted_X->Daughters().front();
    fitted_g2_X    = fitted_X->Daughters().back();

    if(!fitted_Omega || !fitted_g_Omega || !fitted_X || !fitted_g1_X ||!fitted_g1_X ||!fitted_g2_X)
        throw std::runtime_error("Error initializing OmegaEtaG2::MyTreeFitter_t");
}

void OmegaEtaG2::MyTreeFitter_t::HypTestCombis(const TParticleList& photons, doubles& chi2s, doubles& probs, doubles& ggims, doubles& gggims, int& bestIndex)
{

    APLCON::Result_t treefitres;

    bestIndex = -1;
    chi2s = { inf,  inf ,  inf};
    probs = {-inf, -inf , -inf};
    double bestChi2 = inf;

    while(treefitter.NextFit(treefitres)) {

        const auto chi2 = treefitres.Status == APLCON::Result_Status_t::Success ? treefitres.ChiSquare : NaN;
        const auto prob = treefitres.Status == APLCON::Result_Status_t::Success ? treefitres.Probability : NaN;

        const auto combindex = CombIndex(photons, *this);

        chi2s.at(combindex)  = chi2;
        probs.at(combindex)  = prob;
        ggims.at(combindex)  = fitted_X->Get().LVSum.M();
        gggims.at(combindex) = fitted_Omega->Get().LVSum.M();

        if( isfinite(chi2) && chi2 < bestChi2 ) {
            bestIndex = int(combindex);
            bestChi2 = chi2;
        }

    }
}

OmegaEtaG2::MyTreeFitter_t::~MyTreeFitter_t()
{

}


TagChMultiplicity::TagChMultiplicity(HistogramFactory& hf)
{
    const auto tagger = ExpConfig::Setup::GetDetector<TaggerDetector_t>();
    nchannels = tagger->GetNChannels();

    hTagChMult = hf.makeTH1D("Tagger Channel Multiplicity","# hits/event","",BinSettings(10),"tagchmult");


}

void TagChMultiplicity::Fill(const std::vector<TTaggerHit>& t)
{
    vector<double> counts(nchannels, 0);

    for(const auto& thit : t) {
        counts.at(thit.Channel)++;
    }

    for(const auto& m : counts) {
        hTagChMult->Fill(m);
    }
}



void OmegaEtaG2::DebugCounters::TaggerLoopEnd() {
    hPIDLoopsPerTaggerHit->Fill(v.nPIDLoopsPerTaggerHit);
    hKinfitsPerTaggerHit->Fill(v.nKinfitsPerTaggerHit);
    hParticlesFoundPerTH->Fill(v.nParticlesFoundPerTH);
}

void OmegaEtaG2::DebugCounters::AfterKinfitLoop()
{
    ++v.nKinFitsOKPerEvent;
    hHighScoresPerPIDLoop->Fill(v.nHighScoresPerPIDLoop);
    hCoplanarityOKPerPIDLoop->Fill(v.nCoplanarityOKPerPIDLoop);
    hMissingMassOKPerPIDLoop->Fill(v.nMissingMassOKPerPIDLoop);
}

void OmegaEtaG2::DebugCounters::EventEnd()
{
    hTaggerLoops->Fill(v.nTaggerLoops);
    hKinfitsPerEvent->Fill(v.nKinfitsPerEvent);
    hKinFitsOKPerEvent->Fill(v.nKinFitsOKPerEvent);
    hTreeFitsPerEvent->Fill(v.nTreeFitsPerEvent);
    hTaggerHitsAccepted->Fill(v.nTaggerHitAccepted);
}

OmegaEtaG2::DebugCounters::DebugCounters(HistogramFactory& hf)
{
    hTaggerLoops          = hf.makeTH1D("Number of Tagger Loops per Event",       "n Tagger Loops / Event",   "", BinSettings(20), "hTaggerLoops");
    hPIDLoopsPerTaggerHit = hf.makeTH1D("Number of PID Loops per Tagger Hit",     "n PID Loops / Tagger Hit", "", BinSettings(10), "hPIDLoops");
    hKinfitsPerTaggerHit  = hf.makeTH1D("Number of KinFits Loops per Tagger Hit", "n KinFits / Tagger Hit",   "", BinSettings(10), "hKinFitsPerTH");
    hKinfitsPerEvent      = hf.makeTH1D("Number of KinFits Loops per Event",      "n KinFits / Event",        "", BinSettings(10),"hKinFitsPerEvent");
    hKinFitsOKPerEvent    = hf.makeTH1D("Number of KinFits OK per Event",         "n KinFits OK / Event",     "", BinSettings(10), "hKinFitsOKPerEvent");
    hTreeFitsPerEvent     = hf.makeTH1D("Number of TreeFits per Event",           "n Tree Fits / Event",      "", BinSettings(10), "hTreeFitsPerEvent");
    hHighScoresPerPIDLoop = hf.makeTH1D("Number of Kinfit HighScores per PID Loop","n Kin Fit High Scores / PID Loop","", BinSettings(10), "hKinFitHighScoresPerPIDLoop");
    hCoplanarityOKPerPIDLoop = hf.makeTH1D("Coplanarity Cut passed per PID Loop", "n Copl Cut passed / PID Loop","", BinSettings(10), "hCoplanarityOKPerPIDLoop");
    hMissingMassOKPerPIDLoop = hf.makeTH1D("Missing Mass Cut passed per PID Loop", "n MM Cut passed / PID Loop","", BinSettings(10), "hMissingMassOKPerPIDLoop");
    hParticlesFoundPerTH  = hf.makeTH1D("Set of particles found per Tagger Hit",  "n sets / Tagger Hit",        "", BinSettings(10), "hParticlesFoundPerTH");
    hTaggerHitsAccepted  = hf.makeTH1D("Number of Tagger Hits accepted per Event","n Tagger Hits accepted / Event","", BinSettings(10), "hTaggerHitsAccepted");
}


/**
 * @brief Plot Class for the Omega_EtaG analysis
 * @see OmegaEtaG
 */
class OmegaEtaG_Plot : public Plotter {
protected:

    TTree* t = nullptr;

    static const string data_name;

    static bool Contains(const interval<double>& i, const std::vector<double>& d) {
        for(const auto& v : d)
            if(i.Contains(v))
                return true;

        return false;
    }


    static double max(const std::vector<double>& data) {
        return *max_element(data.cbegin(), data.cend());
    }

    static double maxIM(const std::vector<TLorentzVector>& data) {
        return max_element(data.cbegin(), data.cend(), [] (const TLorentzVector& v1, const TLorentzVector& v2) { return v1.M() < v2.M(); })->M();
    }

    static double maxE(const std::vector<TLorentzVector>& data) {
        return max_element(data.cbegin(), data.cend(), [] (const TLorentzVector& v1, const TLorentzVector& v2) { return v1.E() < v2.E(); })->M();
    }

public:

    class OmegaDalitzPlot {

    public:

        static double x(double T1, double T2, double T3) noexcept {
            return (T2 - T1) / (sqrt(2) * (T1+T2+T3));
        }

        static double y(double T1, double T2, double T3) noexcept {
            return T3 / (T1+T2+T3) - 1.0/3.0;
        }

        static vec2 xy(double T1, double T2, double T3) noexcept {
            return vec2(x(T1,T2,T3),y(T1,T2,T3));
        }

        static std::vector<double> getDalitzT(const std::vector<TLorentzVector>& photons, const TLorentzVector& omega) {

            const auto boost = -omega.BoostVector();

            vector<double> T;
            T.reserve(photons.size());
            for( const auto& g : photons) {
                TLorentzVector lv = g;
                lv.Boost(boost);
                T.push_back(lv.E());
            }

            sort(T.begin(),T.end());

            return T;
        }

    protected:
        vec2 calc() {
            return xy(T.at(0), T.at(1), T.at(2));
        }

        std::vector<double> T;

    public:
        vec2 var;

        OmegaDalitzPlot(const std::vector<TLorentzVector>& photons, const TLorentzVector& omega):
            T(getDalitzT(photons,omega)),
            var(calc())
        {}

        bool Next() noexcept {
            const auto r = std::next_permutation(T.begin(), T.end());
            if(r)
                var = calc();
            return r;
        }

    };

    template<typename Hist_t>
    struct MCTrue_Splitter : plot::cuttree::StackedHists_t<Hist_t> {

        // Hist_t should have that type defined
        using Fill_t = typename Hist_t::Fill_t;

        const decltype (physics::OmegaEtaG2::makeChannels()) Channels;

        MCTrue_Splitter(const HistogramFactory& histFac,
                        const plot::cuttree::TreeInfo_t& treeInfo) :
            plot::cuttree::StackedHists_t<Hist_t>(histFac, treeInfo),
            Channels(physics::OmegaEtaG2::makeChannels())
        {
            using plot::histstyle::Mod_t;

            // TODO: derive this from channel map
            this->GetHist(0, data_name, Mod_t::MakeDataPoints(kBlack));
            this->GetHist(1, "Sig",  Mod_t::MakeLine(kRed, 2));
            this->GetHist(2, "Ref",  Mod_t::MakeLine(kGreen, 2));
            // mctrue is never >=3 (and <9) in tree, use this to sum up all MC and all bkg MC
            // see also Fill()
            this->GetHist(3, "Sum_MC", Mod_t::MakeLine(kBlack, 1));
            this->GetHist(4, "Bkg_MC", Mod_t::MakeFill(kGray+1, -1));

        }

        void Fill(const Fill_t& f) {

            const auto mctrue = f.Tree.Channel();

            using plot::histstyle::Mod_t;

            auto get_bkg_name = [] (const unsigned mctrue) {
                const auto entry = physics::OmegaEtaG2::reaction_channels.channels.find(mctrue);

                if(entry!=physics::OmegaEtaG2::reaction_channels.channels.end())
                    return entry->second.name;

                return string("Unknown Decay");
            };

            using plot::histstyle::Mod_t;
            const Hist_t& hist = mctrue<10 ? this->GetHist(mctrue) :
                                             this->GetHist(mctrue,
                                                           get_bkg_name(mctrue),
                                                           Mod_t::MakeLine(plot::histstyle::color_t::GetLight(mctrue-10), 1, kGray+1)
                                                           );


            hist.Fill(f);

            // handle MC_all and MC_bkg
            if(mctrue>0) {
                this->GetHist(3).Fill(f);
                if(mctrue >= 10)
                    this->GetHist(4).Fill(f);
            }
        }
    };

    // define the structs containing the histograms
    // and the cuts. for simple branch variables, that could
    // be combined...

    struct OmegaHist_t {

        static OptionsPtr opts;

        static map<int,double> mc_scale;

        using Tree_t = physics::OmegaEtaG2::OmegaTree_t;

        struct Fill_t {
            const Tree_t& Tree;

            Fill_t(const Tree_t& t) : Tree(t) {}

            double MCScale() const noexcept {
                auto e = mc_scale.find(int(Tree.Channel));
                if(e == mc_scale.end()) {
                    return 1.0;
                }
                return e->second;
            }

            double Weight() const noexcept {
                return Tree.TaggW * MCScale();
            }

            int iBestIndex() const {
                if(Tree.bestHyp == 2) {
                    return Tree.iBestEta;
                } else if(Tree.bestHyp == 1) {
                    return Tree.iBestPi0;
                } else
                    return -1;
            }

            double BestBachelorE() const {
                return iBestIndex() != -1 ? Tree.BachelorE_fitted().at(size_t(iBestIndex())) : NaN;
            }

            int BachelorIndex() const {
                return iBestIndex();
            }

            int nTAPS() const {
                return int(count_if(
                            Tree.photons_detector().begin(), Tree.photons_detector().end(),
                            [] (const int& d) {return d==OmegaBase::DetTAPS;}))
                        + Tree.p_detector == OmegaBase::DetTAPS ? 1 : 0;
            }
        };

        template <typename Hist>
        using fillfunc_t = std::function<void(Hist*, const Fill_t&)>;

        template <typename Hist>
        struct HistFiller_t {
            fillfunc_t<Hist> func;
            Hist* h;
            HistFiller_t(Hist* hist, fillfunc_t<Hist> f): func(f), h(hist) {}
            void Fill(const Fill_t& data) const {
                func(h, data);
            }
        };

        template <typename Hist>
        struct HistMgr : std::list<HistFiller_t<Hist>> {

            using list<HistFiller_t<Hist>>::list;

            void Fill(const Fill_t& data) const {
                for(auto& h : *this) {
                    h.Fill(data);
                }
            }
        };


        HistMgr<TH1D> h1;
        HistMgr<TH2D> h2;

        HistogramFactory HistFac;

        void AddTH1(const string &title, const string &xlabel, const string &ylabel, const BinSettings &bins, const string &name, fillfunc_t<TH1D> f) {
            h1.emplace_back(HistFiller_t<TH1D>(
                                HistFac.makeTH1D(title, xlabel, ylabel, bins, name),f));
        }

        void AddTH2(const string &title, const string &xlabel, const string &ylabel, const BinSettings &xbins, const BinSettings& ybins, const string &name, fillfunc_t<TH2D> f) {
            h2.emplace_back(HistFiller_t<TH2D>(
                                HistFac.makeTH2D(title, xlabel, ylabel, xbins, ybins, name),f));
        }

        OmegaHist_t(const HistogramFactory& hf, plot::cuttree::TreeInfo_t): HistFac(hf) {

            const auto binScale = opts->Get<double>("BinScale", 0.25);

            const auto Bins = [binScale] (const unsigned bins, const double min, const double max) {
                return BinSettings(unsigned(bins*binScale), min, max);
            };

            const BinSettings BachelorEbins = Bins(500, 0, 500);
            const BinSettings ESumbins      = Bins(1600, 0, 1600);
            const BinSettings Ebins         = Bins(1000, 0, 1000);

            const BinSettings probbins = BinSettings(250, 0,   1);

            const BinSettings IMbins        = Bins(1000,   0, 1000);
            const BinSettings gggIMbins     = Bins(60*4, 680.0, 920.0);
            const BinSettings MMbins        = Bins(1000, 400, 1400);

//            const BinSettings MMgggIMbins_X = Bins( 600,   0, 1200);
//            const BinSettings MMgggIMbins_Y = Bins( 750, 500, 2000);

            const BinSettings pThetaBins = Bins( 125,  0,   50);
            const BinSettings  ThetaBins = Bins( 180,  0,   180);
            const BinSettings pEbins     = Bins( 250,  0, 1000);
//            const BinSettings TaggChBins = BinSettings(47);

            const BinSettings TaggTimeBins   = BinSettings(200, -25, 25);
//            const BinSettings CoplBins   = Bins(300, 0, 30.0);

            const BinSettings zVertexBins = Bins(200,-10,10);

            const BinSettings dalitzBins = Bins(200, -0.4, 0.4);
            const BinSettings evtoEbins  = Bins(150,  0, 8);




            // ====== KinFit =======

            //        AddTH1("KinFitChi2",      "#chi^{2}",             "",       Chi2Bins,   "KinFitChi2",
            //               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.KinFitChi2, f.TaggW());
            //        });

            AddTH1("KinFit Probability",      "probability",             "",       probbins,   "KinFitProb",
                   [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.KinFitProb, f.Weight());
                                                 });

            AddTH2("Proton E_{k} Fitted vs Measured", "E_{k} Fitted [MeV]", "E_{k} Measured [MeV]",  pEbins,   pEbins, "p_E_fit_measure",
                   [] (TH2D* h, const Fill_t& f) {
                h->Fill(f.Tree.p_fitted().E() - ParticleTypeDatabase::Proton.Mass(), f.Tree.p().E() - ParticleTypeDatabase::Proton.Mass(), f.Weight());
            });

            //        AddTH2("Proton E_{k} Fitted vs MC True", "E_{k} Fitted [MeV]", "E_{k} MC True [MeV]",  pEbins,   pEbins, "p_E_fit_True",
            //               [] (TH2D* h, const Fill_t& f) {
            //            h->Fill(f.Tree.p_fitted().E() - ParticleTypeDatabase::Proton.Mass(), f.Tree.p.E() - ParticleTypeDatabase::Proton.Mass(), f.TaggW());
            //        });


            // ======= Values after KinFit ======

            AddTH1("3#gamma IM",      "3#gamma IM [MeV]",     "",       gggIMbins,     "ggg_IM",
                   [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.ggg_fitted().M(), f.Weight());
                                                 });

            AddTH1("3#gamma IM unfitted",      "3#gamma IM [MeV]",     "",       gggIMbins,     "ggg_IM_unf",
                   [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.ggg().M(), f.Weight());
                                                 });

            AddTH2("3#gamma IM",      "3#gamma IM [MeV]",     "TaggCH",       gggIMbins, BinSettings(47),   "ggg_IM_taggch",
                   [] (TH2D* h, const Fill_t& f) { h->Fill(f.Tree.ggg_fitted().M(), f.Tree.TaggCh, f.Weight());
                                                 });

            AddTH1("2#gamma sub-IM",  "2#gamma IM [MeV]",     "",       IMbins,     "gg_IM",
                   [] (TH1D* h, const Fill_t& f) {

                for(const auto& v : f.Tree.ggIM_fitted())
                    h->Fill(v.M(), f.Weight());
            });

            AddTH1("Bachelor Photon Energy",  "E [MeV]",     "",       BachelorEbins,     "bachelorE",
                   [] (TH1D* h, const Fill_t& f) {

                for(const auto& v : f.Tree.BachelorE_fitted())
                    h->Fill(v, f.Weight());
            });

            AddTH1("Bachelor Photon Energy|Best Hyp",  "E [MeV]",     "",       BachelorEbins,     "bachelorE_bestHyp",
                   [] (TH1D* h, const Fill_t& f) {

                if( f.iBestIndex() != -1 ) {
                    h->Fill(f.Tree.BachelorE_fitted().at(f.iBestIndex()), f.Weight());
                }
            });

            AddTH2("Bachelor Photon: #theta vs. E", "E [MeV]", "#theta [#circ]",  pEbins,   ThetaBins, "bachelor_theta_E",
                   [] (TH2D* h, const Fill_t& f) {
                if(f.BachelorIndex() < 0)
                    return;
                h->Fill(f.Tree.photons_fitted().at(f.BachelorIndex()).E(), radian_to_degree(f.Tree.photons_fitted().at(f.BachelorIndex()).Theta()), f.Weight());
            });

            AddTH2("Non-Bachelor Photons: #theta vs. E", "E [MeV]", "#theta [#circ]",  pEbins,   ThetaBins, "non_bachelor_theta_E",
                   [] (TH2D* h, const Fill_t& f) {
                const auto bi = f.BachelorIndex();
                if(bi < 0)
                    return;
                for(size_t i=0; i<f.Tree.photons_fitted().size(); ++i) {
                    if(i != size_t(bi))
                        h->Fill(f.Tree.photons_fitted().at(f.BachelorIndex()).E(), radian_to_degree(f.Tree.photons_fitted().at(i).Theta()), f.Weight());
                }

            });

            AddTH1("2#gamma sub-IM|Best Hyp",  "2#gamma IM [MeV]",     "",       IMbins,     "gg_IM_bestHyp",
                   [] (TH1D* h, const Fill_t& f) {

                if( f.iBestIndex() != -1 ) {
                    h->Fill(f.Tree.ggIM().at(f.iBestIndex()).M(), f.Weight());
                }
            });

            AddTH2("Proton #theta vs. E_{k}", "E_{k} [MeV]", "#theta [#circ]",  pEbins,   pThetaBins, "p_theta_E",
                   [] (TH2D* h, const Fill_t& f) {
                h->Fill(f.Tree.p_fitted().E() - ParticleTypeDatabase::Proton.Mass(), radian_to_degree(f.Tree.p_fitted().Theta()), f.Weight());
            });




            AddTH2("Proton #theta vs. E_{k} with Cut", "E_{k} [MeV]","#theta [#circ]",BinSettings(320,0,1600),BinSettings(180,0,180),"p_theta_E_cuts",
                   [] (TH2D* h,const Fill_t& f) {
               h->Fill(f.Tree.p_fitted().E() - ParticleTypeDatabase::Proton.Mass(), radian_to_degree(f.Tree.p_fitted().Theta()),f.Weight());
            });


            //        AddTH2("Missing Mass / 3#gamma IM", "3#gamma IM [MeV]", "MM [MeV]", IMbins,   MMbins,     "mm_gggIM",
            //               [] (TH2D* h, const Fill_t& f) { h->Fill(f.Tree.ggg_fitted().M(), f.Tree.mm().M(), f.TaggW());
            //        });


//            if(opts->Get<bool>("enablePSA",false)) {
//                const auto PSAABins = Bins(  60, 20,   60);
//                const auto PSARBins = Bins( 100,  0,  450);

//                AddTH2("Proton PSA", "PSA Angle [#circ]", "PSA Radius",             PSAABins, PSARBins,   "p_PSA",
//                       [] (TH2D* h, const Fill_t& f) {
//                    h->Fill(f.Tree.p_PSA_Angle, f.Tree.p_PSA_Radius, f.TaggW());
//                });
//            }



    //        AddTH2("3#gamma IM vs 2#gamma IM", "3#gamma IM [MeV]", "max(2#gamma IM) [MeV]", IMbins, IMbins, "ggg_max_gg",
    //               [] (TH2D* h, const Fill_t& f) {
    //            h->Fill(f.Tree.ggg_fitted().M(), maxIM(f.Tree.ggIM_fitted()), f.TaggW());
    //        });

    //        AddTH2("3#gamma E vs 2#gamma IM", "3#gamma E [MeV]", "max(2#gamma IM) [MeV]", IMbins, IMbins, "gggE_max_gg",
    //               [] (TH2D* h, const Fill_t& f) {
    //            h->Fill(f.Tree.ggg_fitted().E()-ParticleTypeDatabase::Omega.Mass(), maxIM(f.Tree.ggIM_fitted()), f.TaggW());
    //        });

    //        AddTH2("3#gamma IM vs 2#gamma max E", "3#gamma IM [MeV]", "max(2#gamma E) [MeV]", IMbins, IMbins, "ggg_max_ggE",
    //               [] (TH2D* h, const Fill_t& f) {
    //            h->Fill(f.Tree.ggg_fitted().M(), maxE(f.Tree.ggIM_fitted()), f.TaggW());
    //        });


            if(opts->Get<bool>("enable-Dalitz",false)) {
                    AddTH2("Dalitz","X","Y", dalitzBins, dalitzBins, "dalitz",
                           [] (TH2D* h, const Fill_t& f) {

                        OmegaDalitzPlot p(f.Tree.photons_fitted(), f.Tree.ggg_fitted());
                        do {
                            h->Fill(p.var.x, p.var.y);
                        } while (p.Next());
                    });
            }


            // ===== Tree Fit =====

            AddTH1("#pi^{0} Hyp, prob", "prob_{#pi^{0}}", "", probbins, "pi0hyp_prob",
                   [] (TH1D* h, const Fill_t& f) {
                const auto& i = f.Tree.iBestPi0;
                if(i >= 0)
                    h->Fill(f.Tree.pi0prob().at(size_t(i)), f.Weight());
            });

            AddTH1("#eta Hyp, prob", "#chi^{2}_{#eta}", "", probbins, "etahyp_prob",
                   [] (TH1D* h, const Fill_t& f) {
                const auto& i = f.Tree.iBestEta;
                if(i >= 0)
                    h->Fill(f.Tree.etaprob().at(size_t(i)), f.Weight());
            });

            AddTH1("#eta Hyp, #omega IM", "m(#omega_{#eta})", "", gggIMbins, "etahyp_omega",
                   [] (TH1D* h, const Fill_t& f) {
                const auto& i = f.Tree.iBestEta;
                if(i >= 0)
                    h->Fill(f.Tree.eta_omega_im().at(size_t(i)), f.Weight());
            });

            AddTH1("#pi^{0} Hyp, #omega IM", "m(#omega_{#pi^{0}}})", "", gggIMbins, "pi0hyp_omega",
                   [] (TH1D* h, const Fill_t& f) {
                const auto& i = f.Tree.iBestPi0;
                if(i >= 0)
                    h->Fill(f.Tree.pi0_omega_im().at(size_t(i)), f.Weight());
            });



            // ====== Crosscheck plots =======

            AddTH2("dEEproton", "E [MeV]", "dE [MeV]", Ebins, evtoEbins, "dEE_proton",
                   [] (TH2D* h, const Fill_t& f) {
                h->Fill(f.Tree.p_fitted().Energy() - ParticleTypeDatabase::Proton.Mass(), f.Tree.p_vetoE);
            });

            AddTH2("dEEphoton", "E [MeV]", "dE [MeV]", Ebins, evtoEbins, "dEE_photon",
                   [] (TH2D* h, const Fill_t& f) {
                for(size_t i=0; i<f.Tree.photons_fitted().size(); ++i) {
                    h->Fill(f.Tree.photons_fitted().at(i).Energy(), f.Tree.photons_vetoE().at(i));
                }
            });


            // ===== Entry Cuts ======

                    AddTH1("Energy of dropped clusters, relative", "dropped E / used E", "", Bins(100,0,.25), "droppedErel",
                            [] (TH1D* h, const Fill_t& f) {

                        h->Fill(f.Tree.CandsunUsedE / f.Tree.CandsUsedE, f.Weight());
                    });

                    AddTH1("Energy of dropped clusters", "dropped E", "", Bins(100,0.0,500.0), "droppedE",
                            [] (TH1D* h, const Fill_t& f) {

                        h->Fill(f.Tree.CandsunUsedE, f.Weight());
                    });

            AddTH1("Missing Mass",      "MM [MeV]",     "",       MMbins,     "mm",
                   [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.mm().M(), f.Weight());
                                                 });

            //        AddTH1("Tagger Channels", "Channel",              "# hits", TaggChBins, "TaggCh",
            //               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.TaggCh, f.TaggW());
            //        });

            //        AddTH1("Coplanarity Angle", "Coplanarity angle [#circ]", "", CoplBins, "CoplAngle",
            //               [] (TH1D* h, const Fill_t& f) { h->Fill(radian_to_degree(f.Tree.copl_angle()), f.TaggW());
            //                                             });

            AddTH1("Tagger Time - CB Average Time", "t [ns]", "",       TaggTimeBins,   "TaggTime",
                   [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.TaggT - f.Tree.CBAvgTime);
                                                 });

            AddTH1("Z Vertex", "z [cm]", "",       zVertexBins,   "zVertex",
                   [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.zVertex, f.Weight());
                                                 });

            AddTH1("Touches Hole Clusters", "n Clusters", "",       BinSettings(5),   "nTouchesHole",
                   [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.nTouchesHole, f.Weight());
                                                 });

            AddTH1("n Candidates", "n Cands", "",       BinSettings(4,4,8),   "nCandidates",
                   [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.nCandsInput, f.Weight());
                                                 });

    //        AddTH2("Touches Hole vs. Kinfit Prob", "KinFit porb", "nClusters Touche Hole", probbins, BinSettings(5), "nTHolesFitProb",
    //               [] (TH2D* h, const Fill_t& f) {
    //            h->Fill(f.Tree.KinFitProb, f.Tree.nTouchesHole, f.TaggW());
    //        });

            AddTH1("CB ESum", "ESum [MeV]", "",       ESumbins,   "CBESum",
                   [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.CBESum, f.Weight());
                                                 });

//            AddTH2("#omega #theta (fitted)", "m(3#gamma) [MeV]", "#theta",       gggIMbins,  BinSettings(18, 0.0, 180.0), "gggim_fitted_theta",
//                   [] (TH2D* h, const Fill_t& f) {

//                const auto& be = f.Tree.TaggE;
//                const auto gp = LorentzVec(vec3(0,0,be),be) + LorentzVec(vec3(0,0,0), ParticleTypeDatabase::Proton.Mass());

//                auto boosted = LorentzVec(f.Tree.ggg_fitted());
//                boosted.Boost(-gp.BoostVector());

//                h->Fill(boosted.M(), radian_to_degree(boosted.Theta()), f.TaggW());

//             });

            if(opts->Get<bool>("ChannelHists", false)) {
                AddTH1("Channels", "", "",       BinSettings(20),   "channels",
                       [] (TH1D* h, const Fill_t& f) {
                    h->Fill(f.Tree.ChannelString().c_str(),1.0);
                });
            }


            // ===== Pulls =====
            if(opts->Get<bool>("enable-Pulls",false)) {

                const auto pullBins = Bins(100,-5,5);

                AddTH1("Pull, Photon CB E", "", "",       pullBins,   "Pull_Photon_CB_E",
                       [] (TH1D* h, const Fill_t& f) {
                    for(size_t i=0; i < 3; ++i) {
                        if(f.Tree.photons_detector().at(i) == 1)
                            h->Fill(f.Tree.photon_E_pulls().at(i),  f.Weight());
                    }
                });

                AddTH1("Pull, Photon CB Theta", "", "",       pullBins,   "Pull_Photon_CB_Theta",
                       [] (TH1D* h, const Fill_t& f) {
                    for(size_t i=0; i < 3; ++i) {
                        if(f.Tree.photons_detector().at(i) == 1)
                            h->Fill(f.Tree.photon_theta_pulls().at(i),  f.Weight());
                    }
                });

                AddTH1("Pull, Photon CB Phi", "", "",       pullBins,   "Pull_Photon_CB_Phi",
                       [] (TH1D* h, const Fill_t& f) {
                    for(size_t i=0; i < 3; ++i) {
                        if(f.Tree.photons_detector().at(i) == 1)
                            h->Fill(f.Tree.photon_phi_pulls().at(i),  f.Weight());
                    }
                });



                AddTH1("Pull, Photon WeightsTAPS E", "", "",       pullBins,   "Pull_Photon_TAPS_E",
                       [] (TH1D* h, const Fill_t& f) {
                    for(size_t i=0; i < 3; ++i) {
                        if(f.Tree.photons_detector().at(i) == 2)
                            h->Fill(f.Tree.photon_E_pulls().at(i),  f.Weight());
                    }
                });

                AddTH1("Pull, Photon TAPS Theta", "", "",       pullBins,   "Pull_Photon_TAPS_Theta",
                       [] (TH1D* h, const Fill_t& f) {
                    for(size_t i=0; i < 3; ++i) {
                        if(f.Tree.photons_detector().at(i) == 2)
                            h->Fill(f.Tree.photon_theta_pulls().at(i),  f.Weight());
                    }
                });

                AddTH1("Pull, Photon TAPS Phi", "", "",       pullBins,   "Pull_Photon_TAPS_Phi",
                       [] (TH1D* h, const Fill_t& f) {
                    for(size_t i=0; i < 3; ++i) {
                        if(f.Tree.photons_detector().at(i) == 2)
                            h->Fill(f.Tree.photon_phi_pulls().at(i),  f.Weight());
                    }
                });

                AddTH1("Pull, Proton CB Theta", "", "",       pullBins,   "Pull_Proton_CB_Theta",
                       [] (TH1D* h, const Fill_t& f) {
                    if(f.Tree.p_detector == 1)
                        h->Fill(f.Tree.p_theta_pull,  f.Weight());
                });
                AddTH1("Pull, Proton CB Phi", "", "",       pullBins,   "Pull_Proton_CB_Phi",
                       [] (TH1D* h, const Fill_t& f) {
                    if(f.Tree.p_detector == 1)
                        h->Fill(f.Tree.p_phi_pull,  f.Weight());
                });


                AddTH1("Pull, Proton TAPS Theta", "", "",       pullBins,   "Pull_Proton_TAPS_Theta",
                       [] (TH1D* h, const Fill_t& f) {
                    if(f.Tree.p_detector == 2)
                        h->Fill(f.Tree.p_theta_pull,  f.Weight());
                });
                AddTH1("Pull, Proton TAPS Phi", "", "",       pullBins,   "Pull_Proton_TAPS_Phi",
                       [] (TH1D* h, const Fill_t& f) {
                    if(f.Tree.p_detector == 2)
                        h->Fill(f.Tree.p_phi_pull,  f.Weight());
                });

                AddTH1("Pull, Bachelor Photon E", "", "",       pullBins,   "Pull_BachelorProton_E",
                       [] (TH1D* h, const Fill_t& f) {
                    if(f.iBestIndex()!= -1)
                        h->Fill(f.Tree.photon_E_pulls().at(size_t(f.iBestIndex())),  f.Weight());
                });

                AddTH1("Bachelor Photon E xcheck", "", "",       Ebins,   "BachelorPhoton_Echeck",
                       [] (TH1D* h, const Fill_t& f) {
                    if(f.BachelorIndex()!= -1) {
                        TVector3 boost = -f.Tree.ggg_fitted().BoostVector();
                        TLorentzVector x = f.Tree.photons_fitted().at(size_t(f.BachelorIndex()));
                        x.Boost(boost);
                        h->Fill(x.E(),  f.Weight());
                    }

                });
            }

        }

        void Fill(const Fill_t& f) const {

            h1.Fill(f);
            h2.Fill(f);

        }

        std::vector<TH1*> GetHists() const {
            vector<TH1*> v;
            v.reserve(h1.size()+h2.size());
            for(auto& e : h1) {
                v.emplace_back(e.h);
            }
            for(auto& e: h2) {
                v.emplace_back(e.h);
            }
            return v;
        }

        static TCutG* makeDalitzCut() {
            TCutG* c = new TCutG("DalitzCut", 3);
            c->SetPoint(0, 0.0,  0.2);
            c->SetPoint(1, -.22, -.11);
            c->SetPoint(2,  .22, -.11);
            return c;
        }

        static TCutG* dalitzCut;


        struct TreeCuts {

            struct pi0Hyp {
                const double cf;
                pi0Hyp(const double cfcut=0.3): cf(cfcut) {}
                bool operator() (const Fill_t& f) const {
                    if(f.Tree.bestHyp != 1)
                        return false;

                    const auto& iBestPi0 = f.Tree.iBestPi0;
                    const auto& pi0prob  = f.Tree.pi0prob()[iBestPi0];

                    return pi0prob > cf;
                }
            };

            struct etaHyp {
                const double cf;
                const double pi0veto_cf;
                etaHyp(const double cfcut=0.3, const double pi0veto=0.01): cf(cfcut), pi0veto_cf(pi0veto) {}
                bool operator() (const Fill_t& f) const {
                    if(f.Tree.bestHyp != 2)
                        return false;

                    const auto& etaprob  = f.Tree.etaprob()[f.Tree.iBestEta];
                    const auto& iBestPi0 = f.Tree.iBestPi0;

                    return etaprob > cf && (iBestPi0==-1 || f.Tree.pi0prob()[iBestPi0] < pi0veto_cf);
                }
            };

            /**
            * @brief Cut on the omega mass in m(3gamma) spectrum
            * @param f
            * @return
            * @note not useful
            */
            static bool wmasscut(const Fill_t& f) noexcept {
                return interval<double>(700,900).Contains(f.Tree.ggg().M());
            }

            /**
             * @brief Cut: Always return true
             * @return true
             */
            static bool dontcare(const Fill_t&) noexcept {
                return true;
            }

            static bool KinFitProb_MM(const Fill_t& f) noexcept {
                return     f.Tree.KinFitProb >  0.01
                        && f.Tree.mm().M() < 1100.0
                        && f.Tree.mm().M() >  780.0;
            }

            /**
             * @brief Cut: require all photons neutral, proton charged
             */
            static bool dEECut(const Fill_t& f) {
                for(const auto& photonVetoE : f.Tree.photons_vetoE()) {
                    if (photonVetoE > .1) return false;
                }
                return f.Tree.p_vetoE > .25;
            }

            struct LineFct {
                double m;
                double b;

                LineFct(const vec2& p1, const vec2& p2) :
                    m((p2.y - p1.y) / (p2.x - p1.x)),
                    b(p1.y - m * p1.x)
                {}

                double operator() (const double& x) const noexcept { return m*x+b; }
            };

            static bool gg_ggg_line_cut(const Fill_t& f) {
                static const LineFct l({667,459}, {902,641});
                const double dist = 30.0;
                const vec2 x = {f.Tree.ggg_fitted().M(), maxIM(f.Tree.ggIM_fitted())};

                return fabs( l(x.x) - x.y ) < dist;
            }

            static bool DalitzCut(const Fill_t& f) {
                OmegaDalitzPlot p(f.Tree.photons_fitted(), f.Tree.ggg_fitted());
                do {
                    if(!dalitzCut->IsInside(p.var.x, p.var.y))
                        return false;
                } while (p.Next());
                return true;
            }

            /**
             * @brief Cut on number of candidates in TAPS
             */
            struct nTAPS {
                const int n;
                bool operator() (const Fill_t& f) const noexcept { return f.nTAPS() <= n; }
                nTAPS(int N): n(N) {}
            };

            struct cosThetaCMOmega : ant::interval<double> {

                using ant::interval<double>::interval;

                bool operator() (const Fill_t& f) const noexcept {
                    LorentzVec omega = f.Tree.ggg_fitted();
                    const auto bt = LorentzVec({0,0,f.Tree.TaggE},f.Tree.TaggE) + LorentzVec({0,0,0}, ParticleTypeDatabase::Proton.Mass());
                    omega.Boost(-bt.BoostVector());
                    return Contains(cos(omega.Theta()));
                }

            };

            /**
             * @brief Cut on number of clusters touching a "hole"
             */
            struct nTouchingHole {
                const int n;
                bool operator() (const Fill_t& f) const noexcept { return f.Tree.nTouchesHole <= n; }
                nTouchingHole(int N): n(N) {}
            };

            struct nCands {
                const unsigned n;
                nCands(unsigned N):n(N){}
                bool operator() (const Fill_t& f) const noexcept { return f.Tree.nCandsInput == n; }
            };

        };

    };

    OmegaEtaG_Plot(const std::string& name, const WrapTFileInput& input, OptionsPtr opts);
    virtual ~OmegaEtaG_Plot();
    long long GetNumEntries() const override { return t->GetEntries(); }
    void ProcessEntry(const long long entry) override;

    plot::cuttree::Tree_t<MCTrue_Splitter<OmegaHist_t>> signal_hists;
    OmegaHist_t::Tree_t tree;



};



OptionsPtr OmegaEtaG_Plot::OmegaHist_t::opts = nullptr;

const string OmegaEtaG_Plot::data_name = "Data";
TCutG* OmegaEtaG_Plot::OmegaHist_t::dalitzCut = OmegaEtaG_Plot::OmegaHist_t::makeDalitzCut();

OmegaEtaG_Plot::OmegaEtaG_Plot(const string &name, const WrapTFileInput &input, OptionsPtr opts):
    Plotter(name, input, opts)
{

    if(!opts->Get<bool>("MCScale",false)) {
        OmegaHist_t::mc_scale.clear();
    }

    const auto tree_name = opts->Get<string>("Tree","OmegaEtaG2/tree");

    if(!input.GetObject(tree_name,t))
        throw Exception("Input TTree " + tree_name + " not found");

    if(!tree.Matches(t))
        throw Exception("Structure of TTree " + tree_name + " does not match WrapTTree");

    tree.LinkBranches(t);


    const auto cuts = [opts] () {
        {

            using plot::cuttree::MultiCut_t;
            using Fill_t = OmegaHist_t::Fill_t;
            using TreeCuts = OmegaHist_t::TreeCuts;

            plot::cuttree::Cuts_t<Fill_t> cuts;

            const auto probCut = opts->Get<double>("ProbCut", 0.01);
            const auto pi0veto_cf = opts->Get<double>("Pi0Veto", 0.01);
            const auto probCutlambda = [probCut] (const Fill_t& f) { return f.Tree.KinFitProb >  probCut; };

            if(opts->Get<bool>("cut-IM-window", false)) {
                cuts.emplace_back(MultiCut_t<Fill_t>{
                                      {"IMwindow", [] (const Fill_t& f) {
                                           auto x= f.Tree.ggg_fitted().M();
                                           return interval<double>(680.0,920.0).Contains(x);
                                       }}}
                                  );
            };

            if(opts->Get<bool>("n4prob",true)) {
                const auto n4probl = [&probCutlambda] (const Fill_t& f) { return probCutlambda(f) && TreeCuts::nCands(4)(f); };
                cuts.emplace_back(MultiCut_t<Fill_t>{
                                      {"n==4+Prob",  n4probl}
                                  });
            } else {
                if(opts->Get<bool>("early4cands", true)) {
                    cuts.emplace_back(MultiCut_t<Fill_t>{
                                          {"4candidates",   TreeCuts::nCands(4) }
                                      });
                }
                LOG(INFO) << "Probability cut: " << probCut;
                cuts.emplace_back(MultiCut_t<Fill_t>{
                                      {"Prob+mm",  [probCut] (const Fill_t& f) { return f.Tree.KinFitProb >  probCut; } }
                                  });

                if(opts->Get<bool>("cut-dEE", false)) {
                    cuts.emplace_back(MultiCut_t<Fill_t>{
                                          {"dEECut",   TreeCuts::dEECut },
                                          {"nodEECut", TreeCuts::dontcare }
                                      });
                }
            }


            if(opts->Get<bool>("enable-cut-TouchesHole", false)) {
                cuts.emplace_back(MultiCut_t<Fill_t>{
                                      {"NoTouchHole",  TreeCuts::nTouchingHole(0)},
                                      {"1TouchHole",   TreeCuts::nTouchingHole(1)},
                                      {"2TouchHole",   TreeCuts::nTouchingHole(2)},
                                      {"3TouchHole",   TreeCuts::nTouchingHole(3)},
                                      {"DontCare",     TreeCuts::dontcare}
                                  });
            }

            if(opts->Get<bool>("cut-Hypotheses", true)) {
                cuts.emplace_back(MultiCut_t<Fill_t>{
                                      {"etaHyp", TreeCuts::etaHyp(probCut,pi0veto_cf)},
                                      {"pi0Hyp", TreeCuts::pi0Hyp(probCut)}
                                  });
            }

            if(opts->Get<bool>("cut-nCands",false)) {
                cuts.emplace_back(MultiCut_t<Fill_t>{
                                      {"n==4", TreeCuts::nCands(4) },
                                      {"n==5", TreeCuts::nCands(5) },
                                      {"n==6", TreeCuts::nCands(6) },
                                      {"nAny", TreeCuts::dontcare },
                                  });
            }

            if(opts->Get<bool>("cut-NTAPS", false)) {
                cuts.emplace_back(MultiCut_t<Fill_t>{
                                      {"nTAPS==0", TreeCuts::nTAPS(0)},
                                      {"nTAPS<=1", TreeCuts::nTAPS(1)},
                                      {"nTAPS<=2", TreeCuts::nTAPS(2)},
                                  });
            }

            if(opts->Get<bool>("cut-cosThetaCM", true)) {
                MultiCut_t<Fill_t> bins;
                const auto nBins = opts->Get<int>("nCosTBins", 40);
                const auto binW  = 2.0/nBins;
                for (int i=0; i<nBins; ++i) {
                    const auto binLowEdge = -1 + i*binW;
                    const auto binHighEdge = binLowEdge + binW;
                    bins.emplace_back(
                                formatter() << "cosT_" << i,
                                TreeCuts::cosThetaCMOmega(binLowEdge, binHighEdge)
                                );
                }
                cuts.emplace_back(std::move(bins));

            }

            if(opts->Get<bool>("cut-PhotonEnergy",false)) {
                cuts.emplace_back(MultiCut_t<Fill_t>{
                                      {"Eg > 50", [] (const Fill_t& f) {
                                           for(const auto& p : f.Tree.photons()) {
                                               if(p.E() < 50.0)
                                                    return false;
                                           }
                                           return true;
                                       }}
                                  });
            }



            if(opts->Get<bool>("cut-ParticleTheta",false)) {
                cuts.emplace_back(MultiCut_t<Fill_t>{
                                      {"Theta", [] (const Fill_t& f) {
                                           for(const auto& p : f.Tree.photons()) {
                                               if(p.Theta() < degree_to_radian(7.0))
                                               return false;
                                           }
                                           if(f.Tree.p().Theta() < degree_to_radian(7.0)) {
                                            return false;
                                           }
                                           return true;
                                       }}
                                  });
            }

            return cuts;
        }

    };


    OmegaHist_t::opts = opts;
    signal_hists = plot::cuttree::Make<MCTrue_Splitter<OmegaHist_t>>(HistFac,cuts());

}

OmegaEtaG_Plot::~OmegaEtaG_Plot() {}

void OmegaEtaG_Plot::ProcessEntry(const long long entry) {
    t->GetEntry(entry);

    plot::cuttree::Fill<MCTrue_Splitter<OmegaHist_t>>(signal_hists, {tree});

}

double OmegaMCCrossSection::EgToW(const double Eg) {
    const auto target = LorentzVec({0,0,0},ParticleTypeDatabase::Proton.Mass());
    const auto beam   = LorentzVec({0,0,Eg},Eg);
    return (beam+target).M();
};

OmegaMCCrossSection::OmegaMCCrossSection(const string &name, OptionsPtr opts):
    Physics(name, opts),
    mcweighting(HistFac, [&opts] () -> const utils::MCWeighting::item_t& {
        if(opts->Get<string>("Meson","") == "omega")
            return utils::MCWeighting::Omega;
        else if(opts->Get<string>("Meson","") == "pi0")
            return utils::MCWeighting::Pi0;
        throw ExceptionOptionNeeded("Wrong meson");
        }()),
    norm( [&opts] () {
    if(opts->Get<string>("Meson","") == "omega")
        return 6.91090694237991715e-02;
    else if(opts->Get<string>("Meson","") == "pi0")
        return 4.90141090503409238e+00;
    throw ExceptionOptionNeeded("Wrong meson");
    }()),
    opt_save_events(opts->Get<bool>("SaveEvents", false)),
    opt_NoWeight(opts->Get<bool>("NoWeight", false)),
    meson([&opts] () -> const ParticleTypeDatabase::Type& {
    if(opts->Get<string>("Meson","") == "omega")
        return ParticleTypeDatabase::Omega;
    else if(opts->Get<string>("Meson","") == "pi0")
        return ParticleTypeDatabase::Pi0;
    throw std::runtime_error("Wrong meson");
    }()),
        rng(0)
    {
    const BinSettings EBins(47,1420,1580);
    const BinSettings ctBins(90,-1,1);
    const BinSettings Ekbins(320,0,1600);
    const BinSettings tbins(180,0,180);

    const BinSettings Wbins(47, EgToW(1420), EgToW(1580));
    counts   = HistFac.makeTH2D("Omega Counts","E_{#gamma} [MeV]","cos(#theta)_{cm}",EBins,ctBins,"countsE");
    counts_w = HistFac.makeTH2D("Omega Counts","W [MeV]","cos(#theta)_{cm}",Wbins,ctBins,"countsW");
    protonET = HistFac.makeTH2D("Proton","E_k [MeV]","#theta [#circ]", Ekbins, tbins,"protonET");
    photonsET  = HistFac.makeTH2D("Photons","E [MeV]","#theta [#circ]", Ekbins, tbins,"photonET");

    cosThetaCMcounts = HistFac.makeTH1D("Event Counts","cos(#theta)_{cm}","counts", BinSettings(40,-1,1),"mesonCounts");
    cosThetaTaggChMCcounts = HistFac.makeTH2D("Event Counts", "cos(#theta)_{cm}","TaggCH",BinSettings(40,-1,1),BinSettings(47),"mesonCounts_taggch");

    omega_Theta          = HistFac.makeTH1D("Omega Theta","cos(#theta)_{cm}","",BinSettings(40,-1,1),"omega_theta");

    proton_Theta_mc      = HistFac.makeTH2D("Proton Theta","cos(#theta)_{cm}","#theta [#circ]",BinSettings(40,-1,1),BinSettings(180,0,180),"proton_theta_mc");
    proton_Phi_mc        = HistFac.makeTH2D("Proton Phi","cos(#theta)_{cm}","#phi [#circ]",BinSettings(40,-1,1),BinSettings(180,-180,180),"proton_phi_mc");
    proton_E_mc          = HistFac.makeTH2D("Proton E","E_k [MeV]","#theta [#circ]",Ekbins,BinSettings(180,0,180),"proton_E_mc" );

    pi0_Theta_mc         = HistFac.makeTH2D("Pi_0 Theta","cos(#theta)_{cm}","#theta [#circ]",BinSettings(40,-1,1),BinSettings(180,0,180),"pi0_Theta");
    pi0_Phi_mc           = HistFac.makeTH2D("Pi_0 Phi","cos(#theta)_{cm}","#phi [#circ]",BinSettings(40,-1,1),BinSettings(180,-180,180),"pi0_Phi");
    pi0_E_mc             = HistFac.makeTH2D("Pi_0 E","E_k [MeV]","#theta [#circ]",Ekbins,BinSettings(180,0,180),"pi0_E");

    gamma_from_pi0_Theta_mc = HistFac.makeTH2D("Gamma Theta","cos(#theta)_{cm}","#theta [#circ]",BinSettings(40,-1,1),BinSettings(180,0,180),"Gamma_from_pi0_Theta");
    gamma_from_pi0_Phi_mc   = HistFac.makeTH2D("Gamma Phi","cos(#theta)_{cm}","#phi [#circ]",BinSettings(40,-1,1),BinSettings(180,-180,180),"Gamma_from_pi0_Phi");
    gamma_from_pi0_E_mc     = HistFac.makeTH2D("Gamma E","E_k [MeV]","#theta [#circ]",Ekbins,BinSettings(180,0,180),"Gamma_from_pi0_E");

    gamma_from_omega_theta_mc = HistFac.makeTH2D("Gamma Theta","cos(#theta)_{cm}","#theta [#circ]",BinSettings(40,-1,1),BinSettings(180,0,180),"Gamma_from_omega_Theta");
    gamma_from_omega_phi_mc   = HistFac.makeTH2D("Gamma Phi","cos(#theta)_{cm}","#phi [#circ]",BinSettings(40,-1,1),BinSettings(180,-180,180),"Gamma_from_omega_Phi");
    gamma_from_omega_E_mc     = HistFac.makeTH2D("Gamma E","E_k [MeV]","#theta [#circ]",Ekbins,BinSettings(180,0,180),"Gamma_from_omega_E");
}

TParticleTree_t getFirst(const ParticleTypeDatabase::Type& t, const TParticleTree_t& tree) {
    auto node = tree->Get();
    if(node->Type() == t) {
        return tree;
    } else {
        for(const auto& d : tree->Daughters()){
            auto r = getFirst(t,d);
            if(r)
                return r;
        }
    }
    return nullptr;
}




void OmegaMCCrossSection::ProcessEvent(const TEvent &event, manager_t &m)
{

    const auto tree = event.MCTrue().ParticleTree;
    if(tree) {

        const bool is_MC = event.Reconstructed().ID.isSet(TID::Flags_t::MC);

        if(is_MC) {
            // until here, no physics cuts were done (THIS IS IMPORTANT)
            // so we can fill this into our mcWeightingEtaPrime instances
            mcweighting.SetParticleTree(tree);
        }

        const auto omega = utils::ParticleTools::FindParticle(meson,tree);
        const auto proton = utils::ParticleTools::FindParticle(ParticleTypeDatabase::Proton,tree);
        const auto photons = utils::ParticleTools::FindParticles(ParticleTypeDatabase::Photon, tree);


        TParticleTree_t pi0_tree   = nullptr;
        TParticleTree_t omega_tree = nullptr;
        TParticlePtr pi0 = nullptr;
        vector<TParticlePtr> gamma_from_pi0;
        vector<TParticlePtr> gamma_from_omega;


        if(omega) {
            const LorentzVec target = {{0,0,0},ParticleTypeDatabase::Proton.Mass()};
            for(const auto& th : event.MCTrue().TaggerHits) {
                const auto beamtarget = target + th.GetPhotonBeam();

                const auto omega_cm = [beamtarget, omega] () {
                    LorentzVec o = *omega;
                    o.Boost(-beamtarget.BoostVector());
                    return o;
                }();

                const auto costheta = cos(omega_cm.Theta());
                const auto w = opt_NoWeight ? 1.0 : mcweighting.GetN(th.PhotonEnergy, costheta);

                if(opt_save_events) {
                    const auto rw = w / norm;
                    if(rw > 1.0 || rng.Uniform(0.0,1.0) <=   rw) {
                        m.SaveEvent();
                    }
                }

                omega_Theta->Fill(cos(omega->Theta()),w);

                //Fitted Ranges: (-0.275 - -0.175) (-0.425 - -0.325)
                if(costheta > -1 && costheta < 1.0){

                    proton_Theta_mc->Fill(costheta,radian_to_degree(proton->Theta()),w);
                    proton_Phi_mc->Fill(cos(proton->Theta()),radian_to_degree(proton->Phi()),w);
                    proton_E_mc->Fill(proton->Ek(),radian_to_degree(proton->Theta()),w);

                    pi0_tree     = getFirst(ParticleTypeDatabase::Pi0,tree);
                    omega_tree   = getFirst(ParticleTypeDatabase::Omega,tree);
                    pi0 = pi0_tree->Get();

                    pi0_Theta_mc->Fill(costheta,radian_to_degree(proton->Theta()),w);
                    pi0_Phi_mc->Fill(cos(pi0->Theta()),radian_to_degree(pi0->Phi()),w);
                    pi0_E_mc->Fill(pi0->Ek(),radian_to_degree(pi0->Theta()),w);

                    if(pi0_tree->Daughters().size() == 2){
                        gamma_from_pi0.push_back(pi0_tree->Daughters().front()->Get());
                        gamma_from_pi0.push_back(pi0_tree->Daughters().back()->Get());
                    }
                    for(const auto& gamma : gamma_from_pi0){
                        gamma_from_pi0_Theta_mc->Fill(costheta,radian_to_degree(gamma->Theta()),w);
                        gamma_from_pi0_Phi_mc->Fill(cos(gamma->Theta()),radian_to_degree(gamma->Phi()),w);
                        gamma_from_pi0_E_mc->Fill(gamma->Ek(),radian_to_degree(gamma->Theta()),w);
                    }


                   if(omega_tree->Daughters().size() == 2){
                       gamma_from_omega.push_back(omega_tree->Daughters().front()->Get());
                       gamma_from_omega.push_back(omega_tree->Daughters().back()->Get());
                   }
                   for(const auto& gamma : gamma_from_omega){
                       if(gamma->M() > 1.0){


                   }
                       else{
                           gamma_from_omega_theta_mc->Fill(costheta,radian_to_degree(gamma->Theta()),w);
                           gamma_from_omega_phi_mc->Fill(cos(gamma->Theta()),radian_to_degree(gamma->Phi()),w);
                           gamma_from_omega_E_mc->Fill(gamma->Ek(),radian_to_degree(gamma->Theta()),w);
                       }

                }
                }


                counts->Fill(th.PhotonEnergy, costheta, w);
                counts_w->Fill(beamtarget.M(), costheta, w);
                protonET->Fill(proton->Ek(), radian_to_degree(proton->Theta()), w);
                for(const auto& p : photons) {
                    photonsET->Fill(p->Ek(), radian_to_degree(p->Theta()), w);
                }
                cosThetaCMcounts->Fill(costheta, w);

                for(const auto& th : event.MCTrue().TaggerHits) {
                    cosThetaTaggChMCcounts->Fill(costheta, th.Channel, w);
                }

                mcweighting.Fill();
            }
        }
    }

}

void OmegaMCCrossSection::ShowResult()
{
    canvas(GetName()) << drawoption("colz") << counts << counts_w <<  protonET << photonsET << proton_Theta_mc << proton_Phi_mc <<proton_E_mc << omega_Theta << pi0_Theta_mc << pi0_Phi_mc << pi0_E_mc << gamma_from_pi0_Theta_mc << gamma_from_pi0_Phi_mc << gamma_from_pi0_E_mc << endc;
}

void OmegaMCCrossSection::Finish()
{
    mcweighting.Finish();
}

map<int,double> OmegaEtaG_Plot::OmegaHist_t::mc_scale = [] () {
    map<int,double> s;
    s[1] = 0.00019;  //sig
    s[2] = 0.084;    //ref
    s[10] = 0.207;   // pi0 -> 2g
    s[13] = 0.418; // eta pi0 -> 4g
    s[11] = 0.418; //2pi0 -> 4g
    s[18] = 0.267; //eta -> 2g
    s[OmegaEtaG2::ReactionChannelList_t::other_index] = 0.3;
    return s;

}();

AUTO_REGISTER_PHYSICS(OmegaMCCrossSection)
AUTO_REGISTER_PHYSICS(OmegaEtaG2)
AUTO_REGISTER_PLOTTER(OmegaEtaG_Plot)
