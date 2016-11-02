#include "TwoPi0_MCSmearing.h"

#include "utils/particle_tools.h"
#include "utils/matcher.h"
#include "utils/Uncertainties.h"
#include "base/Logger.h"
#include "expconfig/ExpConfig.h"
#include "tree/TParticle.h"
#include "tree/TCandidate.h"
#include "tree/TCluster.h"

#include "TH1D.h"
#include "TTree.h"

#include <memory>
#include <cassert>

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace std;



TwoPi0_MCSmearing::TwoPi0_MCSmearing(const string& name, OptionsPtr opts) :
    Physics(name, opts)
{
    auto pi0_range = opts->Get<interval<unsigned>>("nPi0",{1,2});
    if(!pi0_range.IsSane())
        throw runtime_error("Given Pi0 range not sane");

    model = make_shared<utils::UncertaintyModels::FitterSergey>();


    for(unsigned mult=pi0_range.Start();mult<=pi0_range.Stop();mult++) {
        multiPi0.emplace_back(std_ext::make_unique<MultiPi0>(HistFac, mult, model, opts->Get<bool>("SkipFitAndTree", false)));
    }
}

void TwoPi0_MCSmearing::ProcessEvent(const TEvent& event, manager_t&)
{
    const auto& data = event.Reconstructed();
    for(auto& m : multiPi0)
        m->ProcessData(data, event.MCTrue().ParticleTree);
}

void TwoPi0_MCSmearing::ShowResult()
{
    for(auto& m : multiPi0)
        m->ShowResult();
}

void TwoPi0_MCSmearing::Finish()
{
    auto interpolated = dynamic_pointer_cast<const utils::UncertaintyModels::Interpolated>(model);

    if(interpolated) {
        LOG(INFO) << "Interpolated Uncertainty Model Statistics:\n" << *interpolated;
    }
}




TwoPi0_MCSmearing::MultiPi0::MultiPi0(HistogramFactory& histFac, unsigned nPi0, utils::UncertaintyModelPtr FitterModel, bool nofitandnotree) :
    multiplicity(nPi0),
    HistFac(std_ext::formatter() << "m" << multiplicity << "Pi0", histFac, std_ext::formatter() << "m" << multiplicity << "Pi0"),
    nPhotons_expected(multiplicity*2),
    skipfit(nofitandnotree),
    directPi0(getParticleTree(multiplicity)),
    model(FitterModel),
    fitter(std_ext::formatter() << multiplicity << "Pi0", 2*multiplicity, model, true),
    h_missingmass(promptrandom),
    h_fitprobability(promptrandom),
    IM_2g_byMM(promptrandom),
    IM_2g_byFit(promptrandom),
    IM_2g_fitted(promptrandom),
    treefitter("treefit_jusitpi0_"+to_string(nPi0), directPi0, model, true)
{
    fitter.SetZVertexSigma(3.0);
    treefitter.SetZVertexSigma(3.0);


    promptrandom.AddPromptRange({-2.5,2.5});
    promptrandom.AddRandomRange({-50,-5});
    promptrandom.AddRandomRange({  5,50});

    steps = HistFac.makeTH1D("Steps","","#",BinSettings(15),"steps");

    Proton_Coplanarity = HistFac.makeTH1D("p Coplanarity","#delta#phi / degree","",BinSettings(400,-180,180),"Proton_Coplanarity");
    Proton_Angle_True = HistFac.makeTH1D("p Angle to true-matched Rec Proton","#delta#alpha / degree","",BinSettings(400,0,50),"Proton_Angle_True");


    h_missingmass.MakeHistograms(HistFac, "h_missingmass","Missing Mass",BinSettings(400,400, 1400),"MM / MeV","#");
    h_fitprobability.MakeHistograms(HistFac, "fit_probability","KinFitter probability",BinSettings(150,0,1),"p","#");

    BinSettings bins_IM(1400,0,1400);

    IM_2g_byMM.MakeHistograms(HistFac, "IM_2g_byMM","IM 2#gamma by HistFac(multiplicity_str, histFac, multiplicity_str);MM",bins_IM,"IM / MeV","#");
    IM_2g_byFit.MakeHistograms(HistFac, "IM_2g_byFit","IM 2#gamma by Fit",bins_IM,"IM / MeV","#");
    IM_2g_fitted.MakeHistograms(HistFac, "IM_2g_fitted","IM 2#gamma fitted",bins_IM,"IM / MeV","#");

    tree = HistFac.makeTTree("tree");
    t.CreateBranches(tree);

    t.ggIM().resize(nPi0);
    t.photons().resize(nPhotons_expected);
    t.photons_fitted().resize(nPhotons_expected);
    t.photons_PSA().resize(nPhotons_expected);
    t.photons_vetoE().resize(nPhotons_expected);
    t.photons_Time().resize(nPhotons_expected);
    t.fit_photons_E_pulls().resize(nPhotons_expected);
    t.fit_photons_Theta_pulls().resize(nPhotons_expected);
    t.fit_photons_Phi_pulls().resize(nPhotons_expected);

    const auto pion_nodes = treefitter.GetTreeNodes(ParticleTypeDatabase::Pi0);
    assert(pion_nodes.size() == nPi0);

    for(const auto& pion_node : pion_nodes) {
        assert(pion_node->Daughters().size() == 2);
        auto g1 = pion_node->Daughters().front();
        auto g2 = pion_node->Daughters().back();
        pions.emplace_back(make_pair(g1,g2));
    }
}

inline TVector2 getPSAVector(const TParticlePtr& p) {
    if(p->Candidate) {
        const auto cluster = p->Candidate->FindCaloCluster();
        if(cluster) {
            return {cluster->Energy, cluster->ShortEnergy};
        }
    }

    throw std::runtime_error("Incomplete Particle without candiate or CaloCluster");

}

inline unsigned DetectorNum(const Detector_t::Any_t& d) {
    if(d & Detector_t::Type_t::CB) return 1;
    if(d & Detector_t::Type_t::TAPS) return 2;
    return 0;
}

void TwoPi0_MCSmearing::MultiPi0::ProcessData(const TEventData& data, const TParticleTree_t& ptree)
{
    steps->Fill("Seen",1);

    // cut on energy sum and number of candidates

    if(data.Trigger.CBEnergySum <= 550)
        return;
    steps->Fill("CBESum>550MeV",1);

    const auto& cands = data.Candidates;
    const auto nCandidates = cands.size();
    const auto nCandidates_expected = nPhotons_expected+1;
    if(nCandidates != nCandidates_expected)
        return;
    std::string nCandidates_cutstr = std_ext::formatter() << "nCandidates==" << nCandidates_expected;
    steps->Fill(nCandidates_cutstr.c_str(),1);

    // do some MCTrue matching if feasible
    TCandidatePtr proton_mctrue_match = nullptr;
    if(ptree && directPi0 &&
       ptree->IsEqual(directPi0, utils::ParticleTools::MatchByParticleName)) {
        // check if MCTrue matches the found proton
        steps->Fill("Found DirectPi0", 1.0);
        auto true_proton = utils::ParticleTools::FindParticle(ParticleTypeDatabase::Proton, ptree, 1);
        auto true_photons = utils::ParticleTools::FindParticles(ParticleTypeDatabase::Photon, ptree);

        auto mymatcher = [&cands] (const std::vector<TParticlePtr> true_particles) {
            return utils::match1to1(true_particles,
                                    cands.get_ptr_list(),
                                    [] (const TParticlePtr& p1, const TCandidatePtr& p2) {
                return p1->Angle(*p2);
            }, {0.0, std_ext::degree_to_radian(15.0)});
        };

        vector<TParticlePtr> true_all(true_photons);
        true_all.push_back(true_proton);
        const auto match_all = mymatcher(true_all);

        proton_mctrue_match = utils::FindMatched(match_all, true_proton);
    }

    t.isMC      = data.ID.isSet(TID::Flags_t::MC);
    t.CBAvgTime = data.Trigger.CBTiming;

    // iterate over tagger hits

    for(const TTaggerHit& taggerhit : data.TaggerHits) {

        bool proton_found = false;

        t.Tagg_E  = taggerhit.PhotonEnergy;
        t.Tagg_Ch = taggerhit.Channel;
        t.Tagg_W  = promptrandom.FillWeight();

        steps->Fill("Seen taggerhits",1.0);

        promptrandom.SetTaggerHit(taggerhit.Time);
        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        double    best_p_angle         = std_ext::inf;
        double    best_chi2            = std_ext::inf;

        TParticlePtr  selected_proton;
        TParticleList selected_photons;

        std::vector<utils::Fitter::FitParticle> best_fitParticles;

        // use any candidate as proton, and do the analysis (ignore ParticleID stuff)
        for(auto i_proton : cands.get_iter()) {

            const auto proton = std::make_shared<TParticle>(ParticleTypeDatabase::Proton, i_proton);
            std::vector<TParticlePtr> photons;
            for(auto i_photon : cands.get_iter()) {
                if(i_photon == i_proton)
                    continue;
                photons.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::Photon, i_photon));
            }

            assert(photons.size() == nPhotons_expected);

            LorentzVec photon_sum({0,0,0},0);
            for(const auto& p : photons) {
                photon_sum += *p;
            }

            // proton coplanarity

            const double d_phi = std_ext::radian_to_degree(vec2::Phi_mpi_pi(proton->Phi()-photon_sum.Phi() - M_PI ));
            Proton_Coplanarity->Fill(d_phi);

            const interval<double> Proton_Copl_cut(-19, 19);
            if(!Proton_Copl_cut.Contains(d_phi))
                continue;
            const string copl_str = std_ext::formatter() << "Copl p in " << Proton_Copl_cut;
            steps->Fill(copl_str.c_str(),1);

            // simple missing mass cut
            const LorentzVec beam_target = taggerhit.GetPhotonBeam() + LorentzVec({0, 0, 0}, ParticleTypeDatabase::Proton.Mass());
            const LorentzVec missing = beam_target - photon_sum;
            const double missing_mass = missing.M();

            h_missingmass.Fill(missing_mass);
            const interval<double> MM_cut(850, 1000);
            if(MM_cut.Contains(missing_mass)) {
                const string MM_str = std_ext::formatter() << "MM in " << MM_cut;
                steps->Fill(MM_str.c_str(),1.0);
                utils::ParticleTools::FillIMCombinations([this] (double x) {IM_2g_byMM.Fill(x);},  2, photons);
            }

            if(skipfit)
                continue;

            auto angle_p_calcp = std_ext::radian_to_degree(missing.Angle(proton->p));
            if(angle_p_calcp > 20.0)
                continue;
            steps->Fill("p angle < 20.0#circ", 1.0);

            // more sophisticated fitter
            fitter.SetEgammaBeam(taggerhit.PhotonEnergy);
            fitter.SetProton(selected_proton);
            fitter.SetPhotons(selected_photons);
            auto fit_result = fitter.DoFit();

            const auto chi2dof = fit_result.ChiSquare / fit_result.NDoF;

            if(angle_p_calcp < best_p_angle) {

                t.kinfit_chi2dof = chi2dof;
                t.kinfit_prob    = fit_result.Probability;

                if(chi2dof < best_chi2) {
                    best_chi2 = chi2dof;
                    t.ProtonSelectionMatch = 1;
                } else {
                    t.ProtonSelectionMatch = 0;
                }

                selected_proton  = proton;
                selected_photons = photons;

                proton_found = true;

                t.proton_fitted = *fitter.GetFittedProton();

                t.Tagg_E_fitted   = fitter.GetFittedBeamE();
                t.fit_Tagg_E_pull = fitter.GetBeamEPull();

                best_fitParticles = fitter.GetFitParticles();

                t.fit_proton_E_pull     = best_fitParticles.at(0).GetPulls().at(0);
                t.fit_proton_Theta_pull = best_fitParticles.at(0).GetPulls().at(1);
                t.fit_proton_Phi_pull   = best_fitParticles.at(0).GetPulls().at(2);

                const auto photons_fitted = fitter.GetFittedPhotons();

                assert(nPhotons_expected +1  == best_fitParticles.size());
                assert(photons.size()        == nPhotons_expected);
                assert(photons_fitted.size() == nPhotons_expected);

                for(size_t i=0; i< nPhotons_expected; ++i) {

                    t.photons().at(i)        = *selected_photons.at(i);
                    t.photons_fitted().at(i) = *photons_fitted.at(i);
                    t.photons_PSA().at(i)    = getPSAVector(selected_photons.at(i));
                    t.photons_vetoE().at(i)  = selected_photons.at(i)->Candidate->VetoEnergy;
                    t.photons_Time().at(i)   = selected_photons.at(i)->Candidate->Time;

                    t.fit_photons_E_pulls().at(i)     = best_fitParticles.at(i+1).GetPulls().at(0);
                    t.fit_photons_Theta_pulls().at(i) = best_fitParticles.at(i+1).GetPulls().at(1);
                    t.fit_photons_Phi_pulls().at(i)   = best_fitParticles.at(i+1).GetPulls().at(2);
                }

                steps->Fill("P angle OK",1.0);
            }

        } // Loop proton

        if(proton_found) {


            //kinfit

            t.proton        = *selected_proton;

            t.proton_vetoE  = selected_proton->Candidate->VetoEnergy;
            t.proton_Time   = selected_proton->Candidate->Time;

            t.proton_PSA    = getPSAVector(selected_proton);
            t.proton_det    = DetectorNum(selected_proton->Candidate->Detector);

            const auto& vetoCl = selected_proton->Candidate->FindVetoCluster();
            if(vetoCl) {
                t.proton_vetoCh = vetoCl->CentralElement;
            }

            t.ProtonMCTrueMatches = selected_proton->Candidate == proton_mctrue_match;

            tree->Fill();
        } // end KinFit ok

    } // Loop tagger

}

void TwoPi0_MCSmearing::MultiPi0::ShowResult()
{
    // buffer->ResetBranchAddresses();
    canvas(std_ext::formatter() << "TwoPi0_MCSmearing: " << multiplicity << "Pi0")
            << steps
            << Proton_Coplanarity
            << h_missingmass.subtracted
            << IM_2g_byMM.subtracted
            << h_fitprobability.subtracted
            << IM_2g_byFit.subtracted
            << IM_2g_fitted.subtracted
            << Proton_Angle_True
            << endc;
}

ParticleTypeTree TwoPi0_MCSmearing::MultiPi0::getParticleTree(const unsigned nPi0)
{
    if(nPi0==1) {
        return ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0_2g);
    }
    else if(nPi0==2) {
        return ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::TwoPi0_4g);
    }
    else if(nPi0==3) {
        return ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::ThreePi0_6g);
    }

    throw std::runtime_error("Invalid nPi0 specified");
}

AUTO_REGISTER_PHYSICS(TwoPi0_MCSmearing)
