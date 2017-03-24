#include "JustPi0.h"

#include "utils/ParticleTools.h"
#include "utils/Matcher.h"
#include "utils/Uncertainties.h"
#include "base/Logger.h"
#include "expconfig/ExpConfig.h"
#include "tree/TParticle.h"
#include "tree/TCandidate.h"
#include "tree/TCluster.h"
#include "analysis/utils/uncertainties/FitterSergey.h"
#include "analysis/utils/uncertainties/Interpolated.h"

#include "TH1D.h"
#include "TTree.h"

#include <memory>
#include <cassert>

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace std;



JustPi0::JustPi0(const string& name, OptionsPtr opts) :
    Physics(name, opts)
{
    auto pi0_range = opts->Get<interval<unsigned>>("nPi0",{1,2});
    if(!pi0_range.IsSane())
        throw runtime_error("Given Pi0 range not sane");

    model = make_shared<utils::UncertaintyModels::FitterSergey>();


    for(unsigned mult=pi0_range.Start();mult<=pi0_range.Stop();mult++) {
        multiPi0.emplace_back(std_ext::make_unique<MultiPi0>(HistFac, triggersimu, mult, model, opts->Get<bool>("SkipFitAndTree", false)));
    }
}

void JustPi0::ProcessEvent(const TEvent& event, manager_t&)
{
    triggersimu.ProcessEvent(event);
    if(!triggersimu.HasTriggered())
        return;

    const auto& data = event.Reconstructed();
    for(auto& m : multiPi0)
        m->ProcessData(data, event.MCTrue().ParticleTree);
}

void JustPi0::ShowResult()
{
    for(auto& m : multiPi0)
        m->ShowResult();
}

void JustPi0::Finish()
{
    auto interpolated = dynamic_pointer_cast<const utils::UncertaintyModels::Interpolated>(model);

    if(interpolated) {
        LOG(INFO) << "Interpolated Uncertainty Model Statistics:\n" << *interpolated;
    }
}




JustPi0::MultiPi0::MultiPi0(const HistogramFactory &histFac, const utils::TriggerSimulation &triggersimu_,
                            unsigned nPi0, utils::UncertaintyModelPtr FitterModel, bool nofitandnotree) :
    multiplicity(nPi0),
    HistFac(std_ext::formatter() << "m" << multiplicity << "Pi0", histFac, std_ext::formatter() << "m" << multiplicity << "Pi0"),
    triggersimu(triggersimu_),
    nPhotons_expected(multiplicity*2),
    skipfit(nofitandnotree),
    directPi0(getParticleTree(multiplicity)),
    model(FitterModel),
    fitter(model, true),
    h_missingmass(promptrandom),
    h_fitprobability(promptrandom),
    IM_2g_byMM(promptrandom),
    IM_2g_byFit(promptrandom),
    IM_2g_fitted(promptrandom),
    treefitter(directPi0, model, true)
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

TVector2 getPSAVector(const TParticlePtr& p) {
    if(p->Candidate) {
        const auto cluster = p->Candidate->FindCaloCluster();
        if(cluster) {
            return {cluster->Energy, cluster->ShortEnergy};
        }
    }

    throw std::runtime_error("Incomplete Particle without candiate or CaloCluster");

}

unsigned DetectorNum(const Detector_t::Any_t& d) {
    if(d & Detector_t::Type_t::CB) return 1;
    if(d & Detector_t::Type_t::TAPS) return 2;
    return 0;
}

void JustPi0::MultiPi0::ProcessData(const TEventData& data, const TParticleTree_t& ptree)
{
    steps->Fill("Seen",1);

    // cut on number of candidates

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

        t.reaction = 1; // signal
    } else {
        t.reaction = 2; // bkg
    }

    t.isMC      = data.ID.isSet(TID::Flags_t::MC);

    if(!t.isMC) {
        t.reaction = 0; //data
    }
    t.CBAvgTime = triggersimu.GetRefTiming();

    // iterate over tagger hits

    for(const TTaggerHit& taggerhit : data.TaggerHits) {

        bool kinfit_ok = false;

        t.Tagg_E  = taggerhit.PhotonEnergy;
        t.Tagg_Ch = taggerhit.Channel;
        t.Tagg_W  = promptrandom.FillWeight();

        steps->Fill("Seen taggerhits",1.0);

        promptrandom.SetTaggerTime(triggersimu.GetCorrectedTaggerTime(taggerhit));
        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        double    treefit_best_chi2    = 20.0;  // set to min req chi2
        double    kinfit_best_chi2     = 20.0;

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

                if(i_photon->VetoEnergy>0.1)
                    continue;

                photons.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::Photon, i_photon));
            }

            if(photons.size() != nPhotons_expected)
                continue;

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
            if(angle_p_calcp > 15.0)
                continue;
            steps->Fill("p angle < 15.0#circ", 1.0);

            // more sophisticated fitter
            auto fit_result = fitter.DoFit(taggerhit.PhotonEnergy, proton, photons);


            if(fit_result.Status != APLCON::Result_Status_t::Success)
                continue;

            steps->Fill("Fit OK",1.0);
            h_fitprobability.Fill(fit_result.Probability);
            const interval<double> fitprob_cut(0.8, 1);
            if(fitprob_cut.Contains(fit_result.Probability)) {
                const string fitprob_str = std_ext::formatter() << "Fit p in " << fitprob_cut;
                steps->Fill(fitprob_str.c_str(), 1.0);
                utils::ParticleTools::FillIMCombinations([this] (double x) {IM_2g_byFit.Fill(x);},  2, photons);
                utils::ParticleTools::FillIMCombinations([this] (double x) {IM_2g_fitted.Fill(x);},  2, fitter.GetFittedPhotons());
            }

            const auto chi2dof = fit_result.ChiSquare / fit_result.NDoF;

            if(chi2dof < kinfit_best_chi2) {

                t.kinfit_chi2dof = chi2dof;
                t.kinfit_prob    = fit_result.Probability;

                t.proton        = *proton;
                t.proton_fitted = *fitter.GetFittedProton();

                t.Tagg_E_fitted   = fitter.GetFittedBeamE();
                t.fit_Tagg_E_pull = fitter.GetBeamEPull();

                t.proton_vetoE  = proton->Candidate->VetoEnergy;
                t.proton_Time   = proton->Candidate->Time;

                t.proton_PSA    = getPSAVector(proton);
                t.proton_det    = DetectorNum(proton->Candidate->Detector);

                const auto& vetoCl = proton->Candidate->FindVetoCluster();
                if(vetoCl) {
                    t.proton_vetoCh = vetoCl->CentralElement;
                }

                t.ProtonMCTrueMatches = proton->Candidate == proton_mctrue_match;

                best_fitParticles = fitter.GetFitParticles();

                t.fit_proton_E_pull     = best_fitParticles.at(0).GetPulls().at(0);
                t.fit_proton_Theta_pull = best_fitParticles.at(0).GetPulls().at(1);
                t.fit_proton_Phi_pull   = best_fitParticles.at(0).GetPulls().at(2);

                const auto photons_fitted = fitter.GetFittedPhotons();

                assert(nPhotons_expected +1  == best_fitParticles.size());
                assert(photons.size()        == nPhotons_expected);
                assert(photons_fitted.size() == nPhotons_expected);

                for(size_t i=0; i< nPhotons_expected; ++i) {

                    t.photons().at(i)        = *photons.at(i);
                    t.photons_fitted().at(i) = *photons_fitted.at(i);
                    t.photons_PSA().at(i)    = getPSAVector(photons.at(i));
                    t.photons_vetoE().at(i)  = photons.at(i)->Candidate->VetoEnergy;
                    t.photons_Time().at(i)   = photons.at(i)->Candidate->Time;

                    t.fit_photons_E_pulls().at(i)     = best_fitParticles.at(i+1).GetPulls().at(0);
                    t.fit_photons_Theta_pulls().at(i) = best_fitParticles.at(i+1).GetPulls().at(1);
                    t.fit_photons_Phi_pulls().at(i)   = best_fitParticles.at(i+1).GetPulls().at(2);
                }

                selected_proton  = proton;
                selected_photons = photons;

                kinfit_ok = true;

                steps->Fill("Chi2 OK",1.0);
            }

        } // Loop proton

        if(kinfit_ok) {

            t.treefit_prob    = std_ext::NaN;
            t.treefit_chi2dof = std_ext::NaN;

            treefitter.PrepareFits(taggerhit.PhotonEnergy, selected_proton, selected_photons);

            APLCON::Result_t treefitres;
            while(treefitter.NextFit(treefitres)) {

                const double treefit_chi2dof   = treefitres.Status == APLCON::Result_Status_t::Success ? treefitres.ChiSquare   : std_ext::NaN;
                const double treefit_prob      = treefitres.Status == APLCON::Result_Status_t::Success ? treefitres.Probability : std_ext::NaN;

                if(treefit_chi2dof < treefit_best_chi2) {

                    treefit_best_chi2 = treefit_chi2dof;

                    t.treefit_prob    = treefit_prob;
                    t.treefit_chi2dof = treefit_chi2dof;

                    // Fill stuff
                    assert(pions.size() == t.ggIM().size());
                    for(size_t i=0; i<pions.size(); ++i) {
                        LorentzVec pion = *(selected_photons.at(pions.at(i).first->Get().PhotonLeafIndex))
                                          +  *(selected_photons.at(pions.at(i).second->Get().PhotonLeafIndex));
                        t.ggIM().at(i) = pion.M();

                    }

                }
            }

            if(ptree) {
                t.reaction_string = utils::ParticleTools::GetDecayString(ptree);
            } else {
                t.reaction_string = "";
            }

            tree->Fill();
        } // end KinFit ok

    } // Loop tagger

}

void JustPi0::MultiPi0::ShowResult()
{
    // buffer->ResetBranchAddresses();
    canvas(std_ext::formatter() << "JustPi0: " << multiplicity << "Pi0")
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

ParticleTypeTree JustPi0::MultiPi0::getParticleTree(const unsigned nPi0)
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

AUTO_REGISTER_PHYSICS(JustPi0)
