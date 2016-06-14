#include "JustPi0.h"

#include "utils/particle_tools.h"
#include "utils/matcher.h"
#include "utils/FitterUncertainties.h"

#include "expconfig/ExpConfig.h"

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
    for(unsigned mult=1;mult<=opts->Get<unsigned>("nPi0",1);mult++) {
        multiPi0.emplace_back(std_ext::make_unique<MultiPi0>(HistFac, mult, opts->Get<bool>("SkipFitAndTree", false)));
    }
}

void JustPi0::ProcessEvent(const TEvent& event, manager_t&)
{
    const auto& data = event.Reconstructed();
    for(auto& m : multiPi0)
        m->ProcessData(data, event.MCTrue().ParticleTree);
}

void JustPi0::ShowResult()
{
    for(auto& m : multiPi0)
        m->ShowResult();
}




JustPi0::MultiPi0::MultiPi0(HistogramFactory& histFac, unsigned nPi0, bool nofitandnotree) :
    multiplicity(nPi0),
    skipfit(nofitandnotree),
    directPi0(getParticleTree(multiplicity)),
    model(utils::UncertaintyModels::MCExtracted::makeAndLoad()),
    fitter(std_ext::formatter() << multiplicity << "Pi0", 2*multiplicity, model),
    h_missingmass(promptrandom),
    h_fitprobability(promptrandom),
    IM_2g_byMM(promptrandom),
    IM_2g_byFit(promptrandom),
    IM_2g_fitted(promptrandom),
    treefitter("treefit_jusitpi0_"+to_string(nPi0), directPi0, nPi0, model)
{
    std::string multiplicity_str = std_ext::formatter() << "m" << multiplicity << "Pi0";
    HistogramFactory HistFac(multiplicity_str, histFac, multiplicity_str);

    promptrandom.AddPromptRange({-2.5,2.5});
    promptrandom.AddRandomRange({-50,-5});
    promptrandom.AddRandomRange({  5,50});

    steps = HistFac.makeTH1D("Steps","","#",BinSettings(15),"steps");

    Proton_Coplanarity = HistFac.makeTH1D("p Coplanarity","#delta#phi / degree","",BinSettings(400,-180,180),"Proton_Coplanarity");
    Proton_Angle_True = HistFac.makeTH1D("p Angle to true-matched Rec Proton","#delta#alpha / degree","",BinSettings(400,0,50),"Proton_Angle_True");


    h_missingmass.MakeHistograms(HistFac, "h_missingmass","Missing Mass",BinSettings(400,400, 1400),"MM / MeV","#");
    h_fitprobability.MakeHistograms(HistFac, "fit_probability","KinFitter probability",BinSettings(150,0,1),"p","#");

    BinSettings bins_IM(1400,0,1400);

    IM_2g_byMM.MakeHistograms(HistFac, "IM_2g_byMM","IM 2#gamma by MM",bins_IM,"IM / MeV","#");
    IM_2g_byFit.MakeHistograms(HistFac, "IM_2g_byFit","IM 2#gamma by Fit",bins_IM,"IM / MeV","#");
    IM_2g_fitted.MakeHistograms(HistFac, "IM_2g_fitted","IM 2#gamma fitted",bins_IM,"IM / MeV","#");


    tree = HistFac.makeTTree("tree");

    tree->Branch("Tagg_W",          addressof(Tagg_W));
    tree->Branch("Tagg_Ch",         addressof(Tagg_Ch));
    tree->Branch("Tagg_E",          addressof(Tagg_E));

    tree->Branch("proton",          addressof(b_proton));
    tree->Branch("proton_fitted",   addressof(b_proton_fitted));
    tree->Branch("proton_PSA",      addressof(b_proton_PSA));
    tree->Branch("proton_PSA1",     addressof(b_proton_PSA1));
    tree->Branch("proton_PSA2",     addressof(b_proton_PSA2));
    tree->Branch("proton_vetoE",    addressof(b_proton_vetoE));

    tree->Branch("treefit_chi2dof", addressof(treefit_chi2dof));
    tree->Branch("treefit_prob",    addressof(treefit_prob));

    tree->Branch("ProtonMCTrueMatches",addressof(ProtonMCTrueMatches));

    fitter.SetupBranches(tree, "Fit");


    buffer = tree->CloneTree();
    buffer->SetName("buffer");
    buffer->SetDirectory(nullptr);
    tree->CopyAddresses(buffer);
}

void JustPi0::MultiPi0::ProcessData(const TEventData& data, const TParticleTree_t& ptree)
{
    const auto nPhotons_expected = multiplicity*2;

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


    // iterate over tagger hits

    for(const TTaggerHit& taggerhit : data.TaggerHits) {
        steps->Fill("Seen taggerhits",1.0);

        promptrandom.SetTaggerHit(taggerhit.Time);
        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        // clear buffer
        buffer->Reset();
        double    best_chi2     = std_ext::inf;
        Long64_t  best_entry    = -1;
        Long64_t  current_entry = 0;

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


            LorentzVec photon_sum(0,0,0,0);
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
            const LorentzVec beam_target = taggerhit.GetPhotonBeam() + LorentzVec(0, 0, 0, ParticleTypeDatabase::Proton.Mass());
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
            fitter.SetEgammaBeam(taggerhit.PhotonEnergy);
            fitter.SetProton(proton);
            fitter.SetPhotons(photons);
            auto fit_result = fitter.DoFit();



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
            if(chi2dof < 20.0)
                continue;

            steps->Fill("Chi2 OK",1.0);

            Tagg_E  = taggerhit.PhotonEnergy;
            Tagg_Ch = taggerhit.Channel;
            Tagg_W  = promptrandom.FillWeight();

            b_proton  = *proton;
            b_proton_fitted = *fitter.GetFittedProton();

            const double b_proton_shortE = proton->Candidate->FindCaloCluster()->ShortEnergy;
            b_proton_vetoE  = proton->Candidate->VetoEnergy;

            b_proton_PSA    = TVector2(b_proton.Energy()        - ParticleTypeDatabase::Proton.Mass(), b_proton_shortE);
            b_proton_PSA1   = TVector2(b_proton_fitted.Energy() - ParticleTypeDatabase::Proton.Mass(), b_proton_shortE);
            b_proton_PSA2   = TVector2(b_proton_PSA1.X(), b_proton_PSA1.X() / b_proton_PSA.X() * b_proton_shortE);


            ProtonMCTrueMatches = proton->Candidate == proton_mctrue_match;

            treefitter.SetEgammaBeam(taggerhit.PhotonEnergy);
            treefitter.SetProton(proton);
            treefitter.SetPhotons(photons);

            APLCON::Result_t treefitres;
            while(treefitter.NextFit(treefitres)) {

                treefit_chi2dof = treefitres.Status == APLCON::Result_Status_t::Success ? treefitres.ChiSquare : std_ext::NaN;
                treefit_prob    = treefitres.Status == APLCON::Result_Status_t::Success ? treefitres.Probability : std_ext::NaN;

                // Fill stuff


                if(treefit_chi2dof < best_chi2) {
                    best_chi2 = treefit_chi2dof;
                    best_entry = current_entry;
                }

                buffer->Fill();
                current_entry++;

            }

        } // Loop proton

        // save best fit
        if(best_entry != -1) {
            buffer->GetEntry(best_entry);
            tree->Fill();
        }

    } // Loop tagger

}

void JustPi0::MultiPi0::ShowResult()
{
    buffer->ResetBranchAddresses();
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
}

AUTO_REGISTER_PHYSICS(JustPi0)
