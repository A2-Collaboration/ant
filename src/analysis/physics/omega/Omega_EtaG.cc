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


OmegaMCTruePlots::PerChannel_t::PerChannel_t(const string& Title, HistogramFactory& hf):
    title(Title)
{
    proton_E_theta  = hf.makeTH2D(title+" Proton" ,"E [MeV]","#theta [#circ]",BinSettings(800,0,1600),BinSettings(180,0,180), title+"_proton_e_theta");
    photons_E_theta = hf.makeTH2D(title+" Photons","E [MeV]","#theta [#circ]",BinSettings(800,0,1600),BinSettings(180,0,180), title+"_photon_e_theta");
}

void OmegaMCTruePlots::PerChannel_t::Show()
{
    canvas("Omega per Channel: "+title) << drawoption("colz") << proton_E_theta << endc;
}

void OmegaMCTruePlots::PerChannel_t::Fill(const TEventData& d)
{
    auto mctrue_particles = utils::ParticleTypeList::Make(d.ParticleTree);
    const auto& protons = mctrue_particles.Get(ParticleTypeDatabase::Proton);
    if(!protons.empty()) {
        const auto& p = protons.at(0);
        proton_E_theta->Fill(p->Ek(), radian_to_degree(p->Theta()));
    }

    for(const auto& photon : mctrue_particles.Get(ParticleTypeDatabase::Photon)) {
        photons_E_theta->Fill(photon->Ek(),  radian_to_degree(photon->Theta()));
    }
}



OmegaMCTruePlots::OmegaMCTruePlots(const std::string& name, OptionsPtr opts):
    Physics(name, opts), multis(HistFac)
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

    auto mctrue_particles = utils::ParticleTypeList::Make(event.MCTrue().ParticleTree);
    multis.Fill(mctrue_particles.GetAll());
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

    canvas photons("OmegaMCTrue photons E Theta");
    photons << drawoption("colz");

    for(auto& entry : channels) {
        photons << entry.second.photons_E_theta;
    }
    photons << endc;
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
    tree->Branch("p",      &proton_vector);
    tree->Branch("omega",  &omega_vector);
    tree->Branch("gamma1", &gamma1_vector);
    tree->Branch("eta",    &eta_vector);
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

        if(photon.VetoEnergy < 0.5)
            return true;

        if(fabs(vec2::Phi_mpi_pi(proton.Phi - photon.Phi - M_PI)) < degree_to_radian(15.0)) {
            return true;
        }

        return false;
    }

    return true;
}

struct MinTracker {
    double v;
    MinTracker(const double& start = std_ext::inf) : v(start) {}
    bool Track(const double& value) { if(value < v) {v = value; return true;} else return false; }
    double operator()() const { return v; }
};

void OmegaEtaG2::Analyse(const TEventData &data, const TEvent& event, manager_t&) {

    dCounters.EventStart();

    t.Channel = reaction_channels.identify(event.MCTrue().ParticleTree);

    if(t.Channel == ReactionChannelList_t::other_index) {
        if(event.MCTrue().ParticleTree!=nullptr) {
            missed_channels->Fill(utils::ParticleTools::GetDecayString(event.MCTrue().ParticleTree).c_str(), 1.0);
        }
    } else {
        found_channels->Fill(t.Channel);
    }

    TH1* steps = stephists.at(t.Channel);

    steps->Fill("Events seen", 1);


    if(data.Candidates.size() < nCandsMin || data.Candidates.size() > nCandsMax)
        return;

    steps->Fill("NCans OK", 1);

    t.nCandsInput = unsigned(data.Candidates.size());

    TCandidatePtrList cands = copy(data.Candidates);

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

    if(!triggersimu.HasTriggered())
        return;
    steps->Fill("Triggered", 1);

    t.CBAvgTime = triggersimu.GetRefTiming();
    if(!isfinite(t.CBAvgTime))
        return;

    steps->Fill("CB Avg Time OK", 1);



    if(data.TaggerHits.size() > 0)
        steps->Fill("6 has TaggHits", 1);


    t.p_true().SetPxPyPzE(0.0,0.0,0.0,-1.0);
    auto mctrue_particles = utils::ParticleTypeList::Make(event.MCTrue().ParticleTree);
    const auto& mctrue_protons = mctrue_particles.Get(ParticleTypeDatabase::Proton);
    if(mctrue_protons.size() == 1) {
        t.p_true = *(mctrue_protons.front());
    }

    tagChMult.Fill(data.TaggerHits);

    auto dcTaggerHitsAccepted = dTaggerHitsAccepted.getHandle();

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


        MinTracker kinfit_best_chi2(opt_kinfit_chi2cut);

        TParticleList selected_photons;
        TParticleList fitted_photons;
        TParticlePtr  selected_proton;
        TParticlePtr  fitted_proton;
        bool fit_ok = false;

        for(auto it_proton = cands.begin(); it_proton != cands.end(); ++it_proton) {

            dCounters.PIDLoopBegin();

            //proton acceptance checks
            if(!ProtonCheck(*it_proton))
                continue;

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
                continue;

            dCounters.ParticlesFound();


            const TParticle ggg(ParticleTypeDatabase::Omega, LVSum(photons.begin(), photons.end()));

            const auto coplanarity_angle = fabs( vec2::Phi_mpi_pi(proton->Phi() - ggg.Phi() - M_PI));


            if(coplanarity_angle > cut_Copl)
                 continue;

            dCounters.CoplanarityOK();

            const LorentzVec beam_target = TagH.GetPhotonBeam() + target; // make global
            const TParticle missing_vector(ParticleTypeDatabase::Proton, beam_target - ggg);

            if(!cut_missing_mass.Contains(missing_vector.M()))
                continue;

            dCounters.MissingMassOK();

            dCounters.KinfitLoopBegin();
            const auto fitres = fitter.DoFit(TagH.PhotonEnergy, proton, photons);

            if(fitres.Status != APLCON::Result_Status_t::Success)
                continue;

            const auto chi2dof = fitres.ChiSquare / fitres.NDoF;

            if(kinfit_best_chi2.Track(chi2dof)) {

                dCounters.KinfitHighscore();

                fit_ok = true;

                // update Tree branches with new best values

                // Kin Fit
                t.KinFitChi2       = chi2dof;
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
                t.copl_angle = coplanarity_angle;

                const auto gggBoost = -ggg.BoostVector();
                const auto gggBoost_fitted = -ggg_fitted.BoostVector();

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

        dCounters.TaggerLoopEnd();
        tree->Fill();

    } // Tagger Loop

    dCounters.EventEnd();

}

utils::UncertaintyModelPtr OmegaEtaG2::getModel(const string& modelname)
{
    if(modelname == "Oli1") {
        return make_shared<utils::UncertaintyModels::Optimized_Oli1>();
    } else if (modelname == "Interpolated") {
        LOG(INFO) << "Using Interpolated Model";
        return utils::UncertaintyModels::Interpolated::makeAndLoad();

    } else if(modelname == "Sergey") {
        auto& setup = ExpConfig::Setup::Get();

        if( string_starts_with(setup.GetName(), "Setup_2007")) {
            LOG(INFO) << "Using Sergey 2007 Model";
            return make_shared<utils::UncertaintyModels::FitterSergey>(utils::UncertaintyModels::FitterSergey::beamtime_t::Eta_2007);
        }

        LOG(INFO) << "Using Sergey 2014 Model";
        return make_shared<utils::UncertaintyModels::FitterSergey>();

    } else if(modelname == "SergeyProton") {
        LOG(INFO) << "Using Sergey Proton Model";
        auto base = getModel("Sergey");
        return make_shared<utils::UncertaintyModels::MeasuredProton>(base);
    }

    throw std::runtime_error("Invalid model name " + modelname);
}

size_t OmegaEtaG2::CombIndex(const TParticleList& orig, const MyTreeFitter_t& f)
{
    for(size_t i=0; i<orig.size(); ++i) {
        if(orig[i] == f.fitted_g_Omega->Get().Leave->Particle) {
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

    m.channels[10] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0_2g),      "#pi^{0}", kYellow};
    m.channels[11] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::TwoPi0_4g),   "#pi^{0} #pi^{0}", kOrange};
    m.channels[12] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::ThreePi0_6g), "#pi^{0} #pi^{0} #pi^{0}",kGreen-9};
    m.channels[13] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0Eta_4g),   "#pi^{0} #eta",kBlue};
    m.channels[14] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Omega_Pi0PiPPiM_2g),"#omega #rightarrow #pi^{0} #pi^{+} #pi^{-}",kMagenta};
    m.channels[15] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Rho_PiPi),          "#rho #rightarrow #pi^{+} #pi^{-}",kTeal};
    m.channels[16] = {ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0PiPi_2gPiPi),    "#pi^{0} (#rightarrow 2 #gamma) #pi^{+} #pi^{-}", kPink};
    m.channels[m.other_index] = ReactionChannel_t(nullptr, "Others", kCyan);

    return m;
}

OmegaEtaG2::OmegaEtaG2(const std::string& name, OptionsPtr opts):
    OmegaBase(opts->Get<string>("Name", name), opts),
    tree(HistFac.makeTTree("tree")),

    cut_ESum(                     opts->Get<double>(                    "CBESum",               600.0)),
    cut_Copl(    degree_to_radian(opts->Get<double>(                    "CoplAngle",             20.0))),
    photon_E_cb(                  opts->Get<decltype(photon_E_cb)>  (   "PhotonECB",        { 50.0, 1000.0})),
    photon_E_taps(                opts->Get<decltype(photon_E_taps)>(   "PhotonETAPS",      { 50.0, 1400.0})),
    proton_theta(degree_to_radian(opts->Get<decltype(proton_theta)> (   "ProtonThetaRange", { 5.0,   45.0}))),
    cut_missing_mass(             opts->Get<decltype(cut_missing_mass)>("MissingMassWindow", interval<double>::CenterWidth(ParticleTypeDatabase::Proton.Mass(), 450.0))),
    opt_kinfit_chi2cut(           opts->Get<double>(                    "KinFit_Chi2Cut",        15.0)),
    opt_FitZVertex(               opts->Get<bool>(                      "KinFit_FitVertex",     true)),
    opt_strict_Vetos(             opts->Get<bool>(                      "Strict_Vetos",         false)),
    opt_z_sigma(                  opts->Get<double>(                    "ZSigma",               3.0)),

    model(getModel(opts->Get<string>("Model", "SergeyProton"))),
    fitter("OmegaEtaG2", 3, model, opt_FitZVertex),
    fitter_pi0(
        ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Omega_gPi0_3g),
        ParticleTypeDatabase::Pi0,
        model,
        opt_FitZVertex,
        opts->Get<bool>(                      "KinFit_FixOmegaMass",     false)
        ),
    fitter_eta(
        ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Omega_gEta_3g),
        ParticleTypeDatabase::Eta,
        model,
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

    promptrandom.AddPromptRange({- 5,   5});
    promptrandom.AddRandomRange({-50, -10});
    promptrandom.AddRandomRange({ 10,  50});

    t.CreateBranches(tree);

    missed_channels = HistFac.makeTH1D("Unlisted Channels","","Total Events seen",BinSettings(20),"unlistedChannels");
    found_channels  = HistFac.makeTH1D("Listed Channels",  "","Total Events seen",BinSettings(20),"listedChannels");

    for(const auto& c : reaction_channels.channels) {

        stephists[c.first] = HistFac.makeTH1D("Steps: " + c.second.name, "", "", BinSettings(14), "steps_" + to_string(c.first));
        stephists[c.first]->SetLineColor(short(c.second.color));

        if(c.first<20)
            found_channels->GetXaxis()->SetBinLabel(int(c.first+1),c.second.name.c_str());
    }

    LOG(INFO) << "Initialized " << GetName() << ":";
    LOG(INFO) << " CBESum Cut           " << cut_ESum                   << " MeV";
    LOG(INFO) << " Coplanarity Cut      " << radian_to_degree(cut_Copl) << " deg";
    LOG(INFO) << " Min. Photon E (CB)   " << photon_E_cb                << " MeV";
    LOG(INFO) << " Min. Photon E (TAPS) " << photon_E_taps              << " MeV";
    LOG(INFO) << " Proton Theta Angle   " << radian_to_degree(proton_theta) << " deg";
    LOG(INFO) << " Missing Mass Window  " << cut_missing_mass           << " MeV";
    LOG(INFO) << " Max. Chi2/dof KinFit " << opt_kinfit_chi2cut;
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

utils::TreeFitter::nodesetup_t::getter fgetter(const bool& fixmass) {

    if(fixmass)
        return {};

    return [] (const ParticleTypeTree& t) { return utils::TreeFitter::nodesetup_t(1.0, (t->Get() == ParticleTypeDatabase::Omega)); };

}

OmegaEtaG2::MyTreeFitter_t::MyTreeFitter_t(const ParticleTypeTree& ttree, const ParticleTypeDatabase::Type& mesonT, utils::UncertaintyModelPtr model, const bool fix_z_vertex, const bool fix_omega_mass):
    treefitter(
        "treefit_"+mesonT.Name(),
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




OmegaMCTruePlots::CBTAPS_Distribution::CBTAPS_Distribution(HistogramFactory& hf): h("CBTAPS_Distr", hf)
{

}

TH1* OmegaMCTruePlots::CBTAPS_Distribution::addMulti(const unsigned n)
{
    auto hist = h.makeTH1D(formatter() << "Multiplicity " << n, "x CB, n-x TAPS", "", n+1, formatter() << "mult_" << n);
    multiplicity.emplace(pair<unsigned,TH1*>({n, hist}));
    return hist;
}

TH1* OmegaMCTruePlots::CBTAPS_Distribution::getMulti(const unsigned n)
{
    auto entry = multiplicity.find(n);

    if(entry == multiplicity.end()) {
        return addMulti(n);
    } else {
        return entry->second;
    }
}

void OmegaMCTruePlots::CBTAPS_Distribution::Fill(const TParticleList& particles)
{
    auto hist = getMulti(particles.size());
    const auto nCB   = count_if(particles.begin(), particles.end(), [] (const TParticlePtr& p) { return degree_to_radian(interval<double>(20.0, 160.0)).Contains(p->Theta()); });
    const auto nTAPS = count_if(particles.begin(), particles.end(), [] (const TParticlePtr& p) { return degree_to_radian(interval<double>(5.0, 20.0)).Contains(p->Theta()); });
    if(nTAPS+nCB == (int)particles.size())
        hist->Fill(nTAPS);
}

class OmegaEtaG_Plot : public Plotter {
protected:
    TTree* t = nullptr;

    TH1D* h_TaggTime = nullptr;

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

            const unsigned mctrue = unsigned(f.Tree.Channel);

            using plot::histstyle::Mod_t;

            auto get_bkg_name = [] (const unsigned mctrue) {
                const auto entry = physics::OmegaEtaG2::reaction_channels.channels.find(int(mctrue));

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

        using Tree_t = physics::OmegaEtaG2::OmegaTree_t;

        struct Fill_t {
            const Tree_t& Tree;

            Fill_t(const Tree_t& t) : Tree(t) {}

            double TaggW() const {
                return Tree.TaggW;
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

        static BinSettings Bins(const unsigned bins, const double min, const double max) {
            return BinSettings(unsigned(bins), min, max);
        }

        HistMgr<TH1D> h1;
        HistMgr<TH2D> h2;

        const BinSettings Ebins         = Bins(1000, 0, 1000);
        const BinSettings BachelorEbins = Bins(500, 0, 500);
        const BinSettings ESumbins      = Bins(1600, 0, 1600);

        const BinSettings Chi2Bins = BinSettings(250, 0,   25);
        const BinSettings probbins = BinSettings(250, 0,   1);

        const BinSettings IMbins        = Bins(1000,   0, 1000);
        const BinSettings gggIMbins     = Bins( 300, 650,  950);
        const BinSettings MMbins        = Bins(1000, 400, 1400);

        const BinSettings MMgggIMbins_X = Bins( 600,   0, 1200);
        const BinSettings MMgggIMbins_Y = Bins( 750, 500, 2000);

        const BinSettings pThetaBins = Bins( 125,  0,   50);
        const BinSettings pEbins     = Bins( 250,  0, 1000);
        const BinSettings PSAABins   = Bins(  60, 20,   60);
        const BinSettings PSARBins   = Bins( 100,  0,  450);
        const BinSettings TaggChBins = BinSettings(47);

        const BinSettings TaggTimeBins   = BinSettings(200, -25, 25);
        const BinSettings CoplBins   = Bins(300, 0, 30.0);

        const BinSettings zVertexBins = Bins(200,-10,10);

        const BinSettings dalitzBins = Bins(200, -0.4, 0.4);
        const BinSettings evtoEbins  = Bins(150,  0, 8);

        const BinSettings pullBins   = Bins(100,-5,5);

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


            // ====== KinFit =======

            //        AddTH1("KinFitChi2",      "#chi^{2}",             "",       Chi2Bins,   "KinFitChi2",
            //               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.KinFitChi2, f.TaggW());
            //        });

            AddTH1("KinFit Probability",      "probability",             "",       probbins,   "KinFitProb",
                   [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.KinFitProb, f.TaggW());
                                                 });

            AddTH2("Proton E_{k} Fitted vs Measured", "E_{k} Fitted [MeV]", "E_{k} Measured [MeV]",  pEbins,   pEbins, "p_E_fit_measure",
                   [] (TH2D* h, const Fill_t& f) {
                h->Fill(f.Tree.p_fitted().E() - ParticleTypeDatabase::Proton.Mass(), f.Tree.p().E() - ParticleTypeDatabase::Proton.Mass(), f.TaggW());
            });

            //        AddTH2("Proton E_{k} Fitted vs MC True", "E_{k} Fitted [MeV]", "E_{k} MC True [MeV]",  pEbins,   pEbins, "p_E_fit_True",
            //               [] (TH2D* h, const Fill_t& f) {
            //            h->Fill(f.Tree.p_fitted().E() - ParticleTypeDatabase::Proton.Mass(), f.Tree.p.E() - ParticleTypeDatabase::Proton.Mass(), f.TaggW());
            //        });


            // ======= Values after KinFit ======

            AddTH1("3#gamma IM",      "3#gamma IM [MeV]",     "",       gggIMbins,     "ggg_IM",
                   [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.ggg_fitted().M(), f.TaggW());
                                                 });

            AddTH1("3#gamma IM unfitted",      "3#gamma IM [MeV]",     "",       gggIMbins,     "ggg_IM_unf",
                   [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.ggg().M(), f.TaggW());
                                                 });

            AddTH1("2#gamma sub-IM",  "2#gamma IM [MeV]",     "",       IMbins,     "gg_IM",
                   [] (TH1D* h, const Fill_t& f) {

                for(const auto& v : f.Tree.ggIM_fitted())
                    h->Fill(v.M(), f.TaggW());
            });

            AddTH1("Bachelor Photon Energy",  "E [MeV]",     "",       BachelorEbins,     "bachelorE",
                   [] (TH1D* h, const Fill_t& f) {

                for(const auto& v : f.Tree.BachelorE_fitted())
                    h->Fill(v, f.TaggW());
            });

            AddTH1("Bachelor Photon Energy|Best Hyp",  "E [MeV]",     "",       BachelorEbins,     "bachelorE_bestHyp",
                   [] (TH1D* h, const Fill_t& f) {

                if( f.iBestIndex() != -1 ) {
                    h->Fill(f.Tree.BachelorE_fitted().at(f.iBestIndex()), f.TaggW());
                }
            });

            AddTH1("2#gamma sub-IM|Best Hyp",  "2#gamma IM [MeV]",     "",       IMbins,     "gg_IM_bestHyp",
                   [] (TH1D* h, const Fill_t& f) {

                if( f.iBestIndex() != -1 ) {
                    h->Fill(f.Tree.ggIM().at(f.iBestIndex()).M(), f.TaggW());
                }
            });

            AddTH2("Proton #theta vs. E_{k}", "E_{k} [MeV]", "#theta [#circ]",  pEbins,   pThetaBins, "p_theta_E",
                   [] (TH2D* h, const Fill_t& f) {
                h->Fill(f.Tree.p_fitted().E() - ParticleTypeDatabase::Proton.Mass(), radian_to_degree(f.Tree.p_fitted().Theta()), f.TaggW());
            });

            //        AddTH2("Missing Mass / 3#gamma IM", "3#gamma IM [MeV]", "MM [MeV]", IMbins,   MMbins,     "mm_gggIM",
            //               [] (TH2D* h, const Fill_t& f) { h->Fill(f.Tree.ggg_fitted().M(), f.Tree.mm().M(), f.TaggW());
            //        });

            //        AddTH2("Proton PSA", "PSA Angle [#circ]", "PSA Radius",             PSAABins, PSARBins,   "p_PSA",
            //               [] (TH2D* h, const Fill_t& f) {
            //            h->Fill(f.Tree.p_PSA_Angle, f.Tree.p_PSA_Radius, f.TaggW());
            //        });

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


            //        AddTH2("Dalitz","X","Y", dalitzBins, dalitzBins, "dalitz",
            //               [] (TH2D* h, const Fill_t& f) {

            //            OmegaDalitzPlot p(f.Tree.photons_fitted(), f.Tree.ggg_fitted());
            //            do {
            //                h->Fill(p.var.x, p.var.y);
            //            } while (p.Next());
            //        });


            // ===== Tree Fit =====

            AddTH1("#pi^{0} Hyp: prob", "prob_{#pi^{0}}", "", probbins, "pi0hyp_prob",
                   [] (TH1D* h, const Fill_t& f) {
                const auto& i = f.Tree.iBestPi0;
                if(i >= 0)
                    h->Fill(f.Tree.pi0prob().at(size_t(i)), f.TaggW());
            });

            AddTH1("#eta Hyp: prob", "#chi^{2}_{#eta}", "", probbins, "etahyp_prob",
                   [] (TH1D* h, const Fill_t& f) {
                const auto& i = f.Tree.iBestEta;
                if(i >= 0)
                    h->Fill(f.Tree.etaprob().at(size_t(i)), f.TaggW());
            });

            AddTH1("#eta Hyp: #omega IM", "m(#omega_{#eta})", "", gggIMbins, "etahyp_omega",
                   [] (TH1D* h, const Fill_t& f) {
                const auto& i = f.Tree.iBestEta;
                if(i >= 0)
                    h->Fill(f.Tree.eta_omega_im().at(size_t(i)), f.TaggW());
            });

            AddTH1("#pi^{0} Hyp: #omega IM", "m(#omega_{#pi^{0}}})", "", gggIMbins, "pi0hyp_omega",
                   [] (TH1D* h, const Fill_t& f) {
                const auto& i = f.Tree.iBestPi0;
                if(i >= 0)
                    h->Fill(f.Tree.pi0_omega_im().at(size_t(i)), f.TaggW());
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

            //        AddTH1("nCands", "# Candidates", "", BinSettings(4,4,8), "nCands",
            //                [] (TH1D* h, const Fill_t& f) {

            //            h->Fill(f.Tree.nCandsInput, f.TaggW());
            //        });

            //        AddTH1("Energy of dropped clusters, relative", "dropped E / used E", "", Bins(100,0,.25), "droppedErel",
            //                [] (TH1D* h, const Fill_t& f) {

            //            h->Fill(f.Tree.CandsunUsedE / f.Tree.CandsUsedE, f.TaggW());
            //        });

            AddTH1("Missing Mass",      "MM [MeV]",     "",       MMbins,     "mm",
                   [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.mm().M(), f.TaggW());
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
                   [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.zVertex);
                                                 });

            AddTH1("Touches Hole Clusters", "n Clusters", "",       BinSettings(5),   "nTouchesHole",
                   [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.nTouchesHole);
                                                 });
    //        AddTH2("Touches Hole vs. Kinfit Prob", "KinFit porb", "nClusters Touche Hole", probbins, BinSettings(5), "nTHolesFitProb",
    //               [] (TH2D* h, const Fill_t& f) {
    //            h->Fill(f.Tree.KinFitProb, f.Tree.nTouchesHole, f.TaggW());
    //        });

            AddTH1("CB ESum", "ESum [MeV]", "",       ESumbins,   "CBESum",
                   [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.CBESum);
                                                 });

            AddTH2("#omega #theta (fitted)", "m(3#gamma) [MeV]", "#theta",       gggIMbins,  BinSettings(18, 0.0, 180.0), "gggim_fitted_theta",
                   [] (TH2D* h, const Fill_t& f) {

                const auto& be = f.Tree.TaggE;
                const auto gp = LorentzVec(vec3(0,0,be),be) + LorentzVec(vec3(0,0,0), ParticleTypeDatabase::Proton.Mass());

                auto boosted = LorentzVec(f.Tree.ggg_fitted());
                boosted.Boost(-gp.BoostVector());

                h->Fill(boosted.M(), radian_to_degree(boosted.Theta()), f.TaggW());

             });
            // ===== Pulls =====


            //        AddTH1("Pull: Photon CB E", "", "",       pullBins,   "Pull_Photon_CB_E",
            //               [] (TH1D* h, const Fill_t& f) {
            //            for(size_t i=0; i < 3; ++i) {
            //                if(f.Tree.photons_detector().at(i) == 1)
            //                    h->Fill(f.Tree.photon_E_pulls().at(i),  f.TaggW());
            //            }
            //        });

            //        AddTH1("Pull: Photon CB Theta", "", "",       pullBins,   "Pull_Photon_CB_Theta",
            //               [] (TH1D* h, const Fill_t& f) {
            //            for(size_t i=0; i < 3; ++i) {
            //                if(f.Tree.photons_detector().at(i) == 1)
            //                    h->Fill(f.Tree.photon_theta_pulls().at(i),  f.TaggW());
            //            }
            //        });

            //        AddTH1("Pull: Photon CB Phi", "", "",       pullBins,   "Pull_Photon_CB_Phi",
            //               [] (TH1D* h, const Fill_t& f) {
            //            for(size_t i=0; i < 3; ++i) {
            //                if(f.Tree.photons_detector().at(i) == 1)
            //                    h->Fill(f.Tree.photon_phi_pulls().at(i),  f.TaggW());
            //            }
            //        });



            //        AddTH1("Pull: Photon TAPS E", "", "",       pullBins,   "Pull_Photon_TAPS_E",
            //               [] (TH1D* h, const Fill_t& f) {
            //            for(size_t i=0; i < 3; ++i) {
            //                if(f.Tree.photons_detector().at(i) == 2)
            //                    h->Fill(f.Tree.photon_E_pulls().at(i),  f.TaggW());
            //            }
            //        });

            //        AddTH1("Pull: Photon TAPS Theta", "", "",       pullBins,   "Pull_Photon_TAPS_Theta",
            //               [] (TH1D* h, const Fill_t& f) {
            //            for(size_t i=0; i < 3; ++i) {
            //                if(f.Tree.photons_detector().at(i) == 2)
            //                    h->Fill(f.Tree.photon_theta_pulls().at(i),  f.TaggW());
            //            }
            //        });

            //        AddTH1("Pull: Photon TAPS Phi", "", "",       pullBins,   "Pull_Photon_TAPS_Phi",
            //               [] (TH1D* h, const Fill_t& f) {
            //            for(size_t i=0; i < 3; ++i) {
            //                if(f.Tree.photons_detector().at(i) == 2)
            //                    h->Fill(f.Tree.photon_phi_pulls().at(i),  f.TaggW());
            //            }
            //        });

            //        AddTH1("Pull: Proton CB Theta", "", "",       pullBins,   "Pull_Proton_CB_Theta",
            //               [] (TH1D* h, const Fill_t& f) {
            //                if(f.Tree.p_detector == 1)
            //                    h->Fill(f.Tree.p_theta_pull,  f.TaggW());
            //        });
            //        AddTH1("Pull: Proton CB Phi", "", "",       pullBins,   "Pull_Proton_CB_Phi",
            //               [] (TH1D* h, const Fill_t& f) {
            //                if(f.Tree.p_detector == 1)
            //                    h->Fill(f.Tree.p_phi_pull,  f.TaggW());
            //        });


            //        AddTH1("Pull: Proton TAPS Theta", "", "",       pullBins,   "Pull_Proton_TAPS_Theta",
            //               [] (TH1D* h, const Fill_t& f) {
            //                if(f.Tree.p_detector == 2)
            //                    h->Fill(f.Tree.p_theta_pull,  f.TaggW());
            //        });
            //        AddTH1("Pull: Proton TAPS Phi", "", "",       pullBins,   "Pull_Proton_TAPS_Phi",
            //               [] (TH1D* h, const Fill_t& f) {
            //                if(f.Tree.p_detector == 2)
            //                    h->Fill(f.Tree.p_phi_pull,  f.TaggW());
            //        });

            //        AddTH1("Pull: Bachelor Photon E", "", "",       pullBins,   "Pull_BachelorProton_E",
            //               [] (TH1D* h, const Fill_t& f) {
            //                if(f.iBestIndex()!= -1)
            //                    h->Fill(f.Tree.photon_E_pulls().at(f.iBestIndex()),  f.TaggW());
            //        });

            //        AddTH1("Bachelor Photon E xcheck", "", "",       Ebins,   "BachelorPhoton_Echeck",
            //               [] (TH1D* h, const Fill_t& f) {
            //                if(f.BachelorIndex()!= -1) {
            //                    TVector3 boost = -f.Tree.ggg_fitted().BoostVector();
            //                    TLorentzVector x = f.Tree.photons_fitted().at(f.BachelorIndex());
            //                    x.Boost(boost);
            //                    h->Fill(x.E(),  f.TaggW());
            //                }

            //        });

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

            /**
            * @brief omega -> eta gamma hypothesis cut with omega -> pi0 gamma veto
            * @param f
            * @return
            */
            static bool etaHypCut(const Fill_t& f) {

                if(f.Tree.bestHyp != 2)
                    return false;

                const auto& etaprob  = f.Tree.etaprob()[f.Tree.iBestEta];
                const auto& iBestPi0 = f.Tree.iBestPi0;

                return etaprob > 0.3 && (iBestPi0==-1 || f.Tree.pi0prob()[iBestPi0] < 0.01);
            }

            /**
             * @brief omega -> pi0 gamma hypothesis cut
             * @param f
             * @return
             */
            static bool pi0HypCut(const Fill_t& f) {
                if(f.Tree.bestHyp != 1)
                    return false;

                const auto& iBestPi0 = f.Tree.iBestPi0;
                const auto& pi0prob  = f.Tree.pi0prob()[iBestPi0];

                return pi0prob > 0.3;
            }

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
             * @brief Always true
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
        };

        // Sig and Ref channel share some cuts...
        static plot::cuttree::Cuts_t<Fill_t> GetCuts() {

            using plot::cuttree::MultiCut_t;

            plot::cuttree::Cuts_t<Fill_t> cuts;

            cuts.emplace_back(MultiCut_t<Fill_t>{
                                  {"Prob+mm",  TreeCuts::KinFitProb_MM}
                              });

    //                cuts.emplace_back(MultiCut_t<Fill_t>{
    //                                         {"dEECut",   TreeCuts::dEECut },
    //                                         {"nodEECut", TreeCuts::dontcare }
    //                                     });

            cuts.emplace_back(MultiCut_t<Fill_t>{
                                  {"NoTouchHole",  [] (const Fill_t& f){
                                       return f.Tree.nTouchesHole == 0;
                                   }},
                                  {"1TouchHole",  [] (const Fill_t& f){
                                       return f.Tree.nTouchesHole <= 1;
                                   }},
                                  {"2TouchHole",  [] (const Fill_t& f){
                                       return f.Tree.nTouchesHole <= 2;
                                   }},
                                  {"3TouchHole",  [] (const Fill_t& f){
                                       return f.Tree.nTouchesHole <= 3;
                                   }},
                                  {"DontCare",  [] (const Fill_t&){
                                       return true;
                                   }},
                              });

            cuts.emplace_back(MultiCut_t<Fill_t>{
                                  {"etaHyp",               TreeCuts::etaHypCut},
                                  {"pi0Hyp",               TreeCuts::pi0HypCut}
                              });

            return cuts;
        }

    };

    OmegaHist_t::Tree_t tree;

    OmegaEtaG_Plot(const std::string& name, const WrapTFileInput& input, OptionsPtr opts):
        Plotter(name, input, opts) {
        if(!input.GetObject("OmegaEtaG2/tree",t))
            throw Exception("Input TTree not found");

        if(!tree.Matches(t))
            throw Exception("TTree structure does not match WrapTTree");

        tree.LinkBranches(t);

        signal_hists = plot::cuttree::Make<MCTrue_Splitter<OmegaHist_t>>(HistFac);

        h_TaggTime = HistFac.makeTH1D("TaggTime", "t [ns]", "", BinSettings(100,-50,50),"taggtime");

    }

    plot::cuttree::Tree_t<MCTrue_Splitter<OmegaHist_t>> signal_hists;

    virtual ~OmegaEtaG_Plot();

    long long GetNumEntries() const override { return t->GetEntries(); }

    void ProcessEntry(const long long entry) override {
        t->GetEntry(entry);

        plot::cuttree::Fill<MCTrue_Splitter<OmegaHist_t>>(signal_hists, {tree});

        h_TaggTime->Fill(tree.TaggT);
    }

    void ShowResult() override {
        canvas(GetName()) << h_TaggTime << endc;
    }
};

const string OmegaEtaG_Plot::data_name = "Data";
TCutG* OmegaEtaG_Plot::OmegaHist_t::dalitzCut = OmegaEtaG_Plot::OmegaHist_t::makeDalitzCut();

OmegaEtaG_Plot::~OmegaEtaG_Plot() {}

AUTO_REGISTER_PHYSICS(OmegaMCTruePlots)
AUTO_REGISTER_PHYSICS(OmegaMCTree)
AUTO_REGISTER_PHYSICS(OmegaEtaG2)
AUTO_REGISTER_PLOTTER(OmegaEtaG_Plot)
