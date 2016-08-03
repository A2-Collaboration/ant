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
    const auto& protons = d.Particles.Get(ParticleTypeDatabase::Proton);
    if(!protons.empty()) {
        const auto& p = protons.at(0);
        proton_E_theta->Fill(p->Ek(), radian_to_degree(p->Theta()));
    }

    for(const auto& photon : d.Particles.Get(ParticleTypeDatabase::Photon)) {
        photons_E_theta->Fill(photon->Ek(),  radian_to_degree(photon->Theta()));
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

    for(size_t i=0; i<cands.size(); ++i) {

        if(i<nCandsMin)
            t.CandsUsedE   += cands.at(i)->CaloEnergy;
        else
            t.CandsunUsedE += cands.at(i)->CaloEnergy;
    }

    cands.resize(nCandsMin);

    t.CBESum = data.Trigger.CBEnergySum;

    if(t.CBESum < cut_ESum)
        return;

    steps->Fill("CBEsum", 1);

    t.CBAvgTime = event.Reconstructed().Trigger.CBTiming;
    if(!isfinite(t.CBAvgTime))
        return;

    steps->Fill("CB Avg Time OK", 1);



    if(data.TaggerHits.size() > 0)
        steps->Fill("6 has TaggHits", 1);

    tagChMult.Fill(data.TaggerHits);

    auto dcTaggerHitsAccepted = dTaggerHitsAccepted.getHandle();

    for(const TTaggerHit& TagH : data.TaggerHits) {

        dCounters.TaggerLoopBegin();

        promptrandom.SetTaggerHit(TagH.Time - t.CBAvgTime);

        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        dCounters.TaggerHitAccepted();

        dcTaggerHitsAccepted.Count();


        t.TaggW  = promptrandom.FillWeight();
        t.TaggE  = TagH.PhotonEnergy;
        t.TaggCh = TagH.Channel;
        t.TaggT  = TagH.Time;


        double kinfit_best_chi2 = opt_kinfit_chi2cut;

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

                if(PhotonCheck(*it_proton))
                    photons.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::Photon, *it_photon));
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


            fitter.SetEgammaBeam(TagH.PhotonEnergy);
            fitter.SetProton(proton);
            fitter.SetPhotons(photons);

            dCounters.KinfitLoopBegin();
            const auto fitres = fitter.DoFit();

            if(fitres.Status != APLCON::Result_Status_t::Success)
                continue;

            const auto chi2dof = fitres.ChiSquare / fitres.NDoF;

            if(chi2dof < kinfit_best_chi2) {

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

                t.p_theta_pull  = fitparticles.at(0).Theta.Pull;
                t.p_phi_pull    = fitparticles.at(0).Phi.Pull;

                for(size_t i=0; i<nphotons; ++i) {
                    t.photons().at(i)            = *(selected_photons.at(i));
                    t.photons_fitted().at(i)     = *(fitted_photons.at(i));
                    t.photons_Time().at(i)       = selected_photons.at(i)->Candidate->Time;
                    t.photons_vetoE().at(i)      = selected_photons.at(i)->Candidate->VetoEnergy;
                    t.photons_PSA().at(i)        = getPSAVector(selected_photons.at(i));
                    t.photons_detector().at(i)   = getDetectorAsInt(selected_photons.at(i)->Candidate->Detector);
                    t.photon_E_pulls().at(i)     = fitparticles.at(i+1).Ek.Pull;
                    t.photon_theta_pulls().at(i) = fitparticles.at(i+1).Theta.Pull;
                    t.photon_phi_pulls().at(i)   = fitparticles.at(i+1).Phi.Pull;
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
                        t.ggIM().at(combindex) = gg.M();
                        const auto g3_boosted = Boost(*g3, gggBoost);
                        t.BachelorE().at(combindex) = g3_boosted.E;
                    }

                    {
                        const auto& g1 = t.photons_fitted().at(comb[0]);
                        const auto& g2 = t.photons_fitted().at(comb[1]);
                        const auto& g3 = t.photons_fitted().at(comb[2]);

                        const auto gg = g1 + g2;

                        t.ggIM_fitted().at(combindex) = gg.M();

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

            fitter_pi0.treefitter.SetEgammaBeam(TagH.PhotonEnergy);
            fitter_pi0.treefitter.SetPhotons(selected_photons);
            fitter_pi0.treefitter.SetProton(selected_proton);

            fitter_pi0.HypTestCombis(selected_photons,
                                     t.pi0chi2,
                                     t.pi0prob,
                                     t.pi0_im,
                                     t.pi0_omega_im,
                                     t.iBestPi0);



            // Kin fit: test eta hypothesis
            dCounters.TreeFit();

            fitter_eta.treefitter.SetEgammaBeam(TagH.PhotonEnergy);
            fitter_eta.treefitter.SetPhotons(selected_photons);
            fitter_eta.treefitter.SetProton(selected_proton);

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
            if(t.iBestPi0 == -1)
                t.bestHyp = 2; // Eta

            // Eta fits failed, but at least one pi0 fit worked -> pi0 wins
            else if(t.iBestEta == -1) {
                t.bestHyp = 1; // Pi0
            }

            // both hypotheses have valid fit results: select lower chi2
            else
                if(t.pi0chi2().at(t.iBestPi0) < t.etachi2().at(t.iBestEta)) {
                    t.bestHyp = 1;  // PI0
                } else {
                    t.bestHyp = 2;  // ETA
                }

        }

        dCounters.TaggerLoopEnd();
        tree->Fill();

        // TODO

        {
//        TParticleList rec_photons(3);
//        TParticlePtr  rec_proton = nullptr;
//        TParticleList true_particles(4);

//        t.ggIM_real = NaN;
//        t.ggIM_comb = {NaN, NaN};

//        const auto& particletree = event.MCTrue().ParticleTree;

//        if(particletree && (t.Channel == 1 || t.Channel == 2)) {
//            particletree->Map_level([&true_particles] (const TParticlePtr& p, const size_t& level) {

//                if(level == 1) {
//                    if(p->Type() == ParticleTypeDatabase::Proton) {
//                        FASSERT(true_particles[3] == nullptr);
//                        true_particles[3] = p;
//                    }
//                }

//                if(p->Type() == ParticleTypeDatabase::Photon) {

//                    if(level==2) {
//                        FASSERT(true_particles[0]==nullptr);
//                        true_particles[0] = p;

//                    } else if(level == 3) {

//                        if(!true_particles[1]) {
//                            true_particles[1] = p;
//                        } else {
//                            FASSERT(true_particles[2]==nullptr);
//                            true_particles[2] = p;
//                        }
//                    }
//                }
//            });

//            FASSERT(true_particles[0]!=nullptr);
//            FASSERT(true_particles[1]!=nullptr);
//            FASSERT(true_particles[2]!=nullptr);
//            FASSERT(true_particles[3]!=nullptr);

//            t.p_true = *true_particles[3];

//            const auto matched  = utils::match1to1(true_particles, data.Particles.GetAll(), TParticle::CalcAngle, {0.0, degree_to_radian(15.0)});

//            if(matched.size() == true_particles.size()) {

//                rec_photons[0] = utils::FindMatched(matched, true_particles[0]);
//                rec_photons[1] = utils::FindMatched(matched, true_particles[1]);
//                rec_photons[2] = utils::FindMatched(matched, true_particles[2]);
//                rec_proton     = utils::FindMatched(matched, true_particles[3]);

//                FASSERT(rec_photons[0]!=nullptr);
//                FASSERT(rec_photons[1]!=nullptr);
//                FASSERT(rec_photons[2]!=nullptr);
//                FASSERT(rec_proton    !=nullptr);


//                t.p_matched = (rec_proton->Type() == ParticleTypeDatabase::Proton);

//                t.ggIM_real      = IM(rec_photons[1], rec_photons[2]);
//                t.ggIM_comb()[0] = IM(rec_photons[0], rec_photons[1]);
//                t.ggIM_comb()[1] = IM(rec_photons[0], rec_photons[2]);

//            }

//        }
        }

    } // Tagger Loop

    dCounters.EventEnd();

}

utils::UncertaintyModelPtr OmegaEtaG2::getModel(const string& modelname)
{
    if(modelname == "Oli1") {
        return make_shared<utils::UncertaintyModels::Optimized_Oli1>();
    } else if (modelname == "Interpolated") {
        return utils::UncertaintyModels::Interpolated::makeAndLoad(
                      make_shared<utils::UncertaintyModels::Optimized_Oli1>(),
                      utils::UncertaintyModels::Interpolated::Mode_t::Fit,
                      false);

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
    m.channels[m.other_index] = ReactionChannel_t(nullptr, "Others", kCyan);

    return m;
}

OmegaEtaG2::OmegaEtaG2(const std::string& name, OptionsPtr opts):
    OmegaBase(opts->Get<string>("Name", name), opts),
    tree(HistFac.makeTTree("tree")),

    cut_ESum(                     opts->Get<double>(                    "CBESum",               550.0)),
    cut_Copl(    degree_to_radian(opts->Get<double>(                    "CoplAngle",             15.0))),
    photon_E_cb(                  opts->Get<decltype(photon_E_cb)>  (   "PhotonECB",        { 50.0, 1600.0})),
    photon_E_taps(                opts->Get<decltype(photon_E_taps)>(   "PhotonETAPS",      {100.0, 1600.0})),
    proton_theta(degree_to_radian(opts->Get<decltype(proton_theta)> (   "ProtonThetaRange", {  4.0,   45.0}))),
    cut_missing_mass(             opts->Get<decltype(cut_missing_mass)>("MissingMassWindow",{780.0, 1200.0})),
    opt_kinfit_chi2cut(           opts->Get<double>(                    "KinFit_Chi2Cut",        10.0)),
    opt_FitZVertex(               opts->Get<bool>(                      "KinFit_FitVertex",     false)),

    model(getModel(opts->Get<string>("Model", "Oli1"))),
    fitter("OmegaEtaG2", 3, model, opt_FitZVertex),
    fitter_pi0(
        ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Omega_gPi0_3g),
        ParticleTypeDatabase::Pi0,
        model,
        opt_FitZVertex
        ),
    fitter_eta(
        ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Omega_gEta_3g),
        ParticleTypeDatabase::Eta,
        model,
        opt_FitZVertex
        ),
    tagChMult(HistFac),
    dTaggerHitsAccepted(HistFac.makeTH1D("Tagger Hits Accepted Per Event","","",BinSettings(10),"dTHAcceptedperEvent")),
    dCounters(HistFac)

{

    promptrandom.AddPromptRange({-5,5});
    promptrandom.AddRandomRange({-20, -10});
    promptrandom.AddRandomRange({ 10,  20});

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
    LOG(INFO) << " Proton Theta Angle   " << proton_theta               << " deg";
    LOG(INFO) << " Missing Mass Window  " << cut_missing_mass           << " MeV";
    LOG(INFO) << " Max. Chi2/dof KinFit " << opt_kinfit_chi2cut;
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



OmegaEtaG2::MyTreeFitter_t::MyTreeFitter_t(const ParticleTypeTree& ttree, const ParticleTypeDatabase::Type& mesonT, utils::UncertaintyModelPtr model, const bool fix_z_vertex):
    treefitter(
        "treefit_"+mesonT.Name(),
        ttree,
        model, fix_z_vertex,
        [] (const ParticleTypeTree& t) { return utils::TreeFitter::nodesetup_t(1.0, (t->Get() == ParticleTypeDatabase::Omega)); }
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


AUTO_REGISTER_PHYSICS(OmegaMCTruePlots)
AUTO_REGISTER_PHYSICS(OmegaMCTree)
AUTO_REGISTER_PHYSICS(OmegaEtaG2)
