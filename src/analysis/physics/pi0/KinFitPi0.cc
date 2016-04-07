#include "KinFitPi0.h"

#include "utils/particle_tools.h"
#include "utils/matcher.h"
#include "base/Logger.h"

#include "expconfig/ExpConfig.h"
#include "base/interval_algo.h"

#include "TH1D.h"
#include "TTree.h"

#include <string>
#include <memory>
#include <cassert>

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace std;



std::shared_ptr<utils::Fitter::UncertaintyModel> KinFitPi0::getModel(const string& model_name) const
{
    if(model_name == "ConstantRelativeE") {
        return utils::UncertaintyModels::ConstantRelativeE::makeMCLongTarget();
    }

    if(model_name == "MCExtracted"){
        return utils::UncertaintyModels::MCExtracted::makeAndLoad();
    }

    if(model_name == "ConstantRelativeEpow") {
        auto s = utils::UncertaintyModels::ConstantRelativeEpow::makeMCLongTarget();
//        s->Eexp_cb   = opts.Get("Eexp_cb",   s->Eexp_cb);
//        s->Eexp_taps = opts.Get("Eexp_taps", s->Eexp_taps);
        return s;
    }

    if(model_name == "Zero") {
        auto s = utils::UncertaintyModels::Constant::make();
        s->photon_cb   = {0,0,0};
        s->photon_taps = {0,0,0};
        s->proton_cb   = {0,0,0};
        s->proton_taps = {0,0,0};
        return s;
    }

    if(model_name == "Theoretical") {
        return std::make_shared<utils::UncertaintyModels::Optimized_Oli1>();
    }

    // fallback
    throw std::runtime_error("no fitter uncertainty model specified");
}

KinFitPi0::KinFitPi0(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    fitter_model(getModel(opts->Get<string>("FitModel", "ConstantRelativeE"))),
    smear_model(getModel(opts->Get<string>("SmearModel", "Zero"))),
    smear(smear_model),
    opt_sweep_param(opts->Get<bool>("Sweep", false))
{

    for(unsigned mult=1;mult<=opts->Get<unsigned>("nPi0",3);mult++) {

        auto mpi0 = std_ext::make_unique<MultiPi0>(HistFac, mult, fitter_model);

        // if using the theoretical model, add some branches to the oputput trees
        if(auto theom = std::dynamic_pointer_cast<utils::UncertaintyModels::Optimized>(fitter_model)) {
            mpi0->tree->Branch("T_cb_photon_theta_const", addressof(theom->cb_photon_theta_const));
            mpi0->tree->Branch("T_cb_photon_theta_Sin",   addressof(theom->cb_photon_theta_Sin));
            mpi0->tree->Branch("T_cb_photon_E_rel",       addressof(theom->cb_photon_E_rel));
            mpi0->tree->Branch("T_cb_photon_E_exp",       addressof(theom->cb_photon_E_exp));
        }

        multiPi0.emplace_back(move(mpi0));

    }

}

void KinFitPi0::ProcessEvent(const TEvent& event, manager_t&)
{

    if(opt_sweep_param) {
        if(auto theom = std::dynamic_pointer_cast<utils::UncertaintyModels::Optimized>(fitter_model)) {

            const interval<double> theta_const_i = { .5, 3.5 };
            const interval<double> theta_Sin_i   = { 3, 8 };
            const unsigned steps = 3;

            for(const auto thsin : Range(theta_Sin_i, steps)) {
                theom->cb_photon_theta_Sin = thsin;

                for(const auto thconst : Range(theta_const_i, steps)) {
                    theom->cb_photon_theta_const = thconst;

                    // fit
                    const auto& data = event.Reconstructed();
                    for(auto& m : multiPi0)
                        m->ProcessData(data, smear);
                }

            }

        }
    } else {
        const auto& data = event.Reconstructed();
        for(auto& m : multiPi0)
            m->ProcessData(data, smear);
    }
}

void KinFitPi0::ShowResult()
{
    for(auto& m : multiPi0)
        m->ShowResult();
}




KinFitPi0::MultiPi0::MultiPi0(HistogramFactory& histFac, unsigned nPi0, std::shared_ptr<const utils::Fitter::UncertaintyModel> model) :
    multiplicity(nPi0),
    HistFac(std_ext::formatter() << "m" << multiplicity << "Pi0", histFac,
            std_ext::formatter() << "m" << multiplicity << "Pi0"),
    IM_perms(BuildIMPerms(multiplicity)),
    fitter(std_ext::formatter() << multiplicity << "Pi0", 2*multiplicity, model),
    h_missingmass(promptrandom),
    h_fitprobability(promptrandom),
    IM_2g_byFit(promptrandom),
    IM_2g_fitted(promptrandom)
{
    promptrandom.AddPromptRange({-2.5,1.5}); // slight offset due to CBAvgTime reference
    promptrandom.AddRandomRange({-50,-10});  // just ensure to be way off prompt peak
    promptrandom.AddRandomRange({  10,50});

    steps = HistFac.makeTH1D("Steps","","#",BinSettings(15),"steps");

    Proton_Coplanarity = HistFac.makeTH1D("p Coplanarity","#delta#phi / degree","",BinSettings(400,-180,180),"Proton_Coplanarity");

    h_taggtime = HistFac.makeTH1D("Tagged Time","t / ns", "", BinSettings(300,-60,60), "h_taggtime");
    h_missingmass.MakeHistograms(HistFac, "h_missingmass","Missing Mass",BinSettings(400,400, 1400),"MM / MeV","#");
    h_fitprobability.MakeHistograms(HistFac, "fit_probability","KinFitter probability",BinSettings(150,0,1),"p","#");

    tree = HistFac.makeTTree("tree");

    BinSettings bins_IM(1400,0,1400);

    IM_2g_byFit.MakeHistograms(HistFac, "IM_2g_byFit","IM 2#gamma by Fit",bins_IM,"IM / MeV","#");
    IM_2g_fitted.MakeHistograms(HistFac, "IM_2g_fitted","IM 2#gamma fitted",bins_IM,"IM / MeV","#");

    fitter.SetupBranches(tree);
    tree->Branch("g1",    &b_g1);
    tree->Branch("g2",    &b_g2);
    tree->Branch("p",     &b_p);
    tree->Branch("TaggW", &b_tagw);

}

void KinFitPi0::MultiPi0::ProcessData(const TEventData& data, const utils::MCSmear& smear)
{
    const auto nPhotons_expected = multiplicity*2;

    steps->Fill("Seen",1);

    const auto& cands = data.Candidates;
    const auto nCandidates = cands.size();
    const auto nCandidates_expected = nPhotons_expected+1;
    if(nCandidates != nCandidates_expected)
        return;
    std::string nCandidates_cutstr = std_ext::formatter() << "nCandidates==" << nCandidates_expected;
    steps->Fill(nCandidates_cutstr.c_str(),1);

    // use any candidate as proton, and do the analysis (ignore ParticleID stuff)

    for(auto i_proton : cands.get_iter()) {

        const auto proton = smear.Smear(std::make_shared<TParticle>(ParticleTypeDatabase::Proton, i_proton));
        std::vector<TParticlePtr> photons;
        unsigned nPhotons_CB = 0;
        unsigned nPhotons_TAPS = 0;
        bool touchesHole = false;
        for(auto i_photon : cands.get_iter()) {
            if(i_photon == i_proton)
                continue;
            photons.emplace_back(smear.Smear(make_shared<TParticle>(ParticleTypeDatabase::Photon, i_photon)));

            if(i_photon->Detector & Detector_t::Type_t::CB)
                nPhotons_CB++;
            else if(i_photon->Detector & Detector_t::Type_t::TAPS)
                nPhotons_TAPS++;

            if(i_photon->FindCaloCluster()->HasFlag(TCluster::Flags_t::TouchesHole))
                touchesHole = true;
        }
        assert(photons.size() == nPhotons_expected);

        if(nPhotons_CB != nPhotons_expected) {
            steps->Fill("Photons in TAPS", 1.0);
            continue;
        }

        if(touchesHole) {
            steps->Fill("TouchesHole", 1.0);
            continue;
        }

        // check invariant mass
        bool im_ok = false;
        for(auto& idx : IM_perms) {
            bool all_inside_pi0 = true;
            for(unsigned i=0;i<multiplicity;i++) {
                const auto im = (*photons[idx[2*i]] + *photons[idx[2*i+1]]).M();
                const auto Pi0 = ParticleTypeDatabase::Pi0.GetWindow(20);
                if(!Pi0.Contains(im)) {
                    all_inside_pi0 = false;
                    break;
                }
            }
            if(all_inside_pi0) {
                im_ok = true;
                break;
            }
        }
        if(!im_ok) {
            steps->Fill("IM not ok",1.0);
            continue;
        }


        LorentzVec photon_sum(0,0,0,0);
        for(const auto& p : photons) {
            photon_sum += *p;
        }

        // proton coplanarity check

        const double d_phi = std_ext::radian_to_degree(vec2::Phi_mpi_pi(proton->Phi()-photon_sum.Phi() - M_PI ));
        Proton_Coplanarity->Fill(d_phi);

        // iterate over tagger hits

        b_g1 = *photons.at(0);
        b_g2 = *photons.at(1);
        b_p  = *proton;

        for(const TTaggerHit& taggerhit : data.TaggerHits) {


            steps->Fill("Seen taggerhits",1.0);
            const auto taggtime = taggerhit.Time - data.Trigger.CBTiming;
            promptrandom.SetTaggerHit(taggtime);
            if(promptrandom.State() == PromptRandom::Case::Outside)
                continue;
            h_taggtime->Fill(taggtime);

            b_tagw = promptrandom.FillWeight();

            // simple missing mass cut
            const LorentzVec beam_target = taggerhit.GetPhotonBeam() + LorentzVec(0, 0, 0, ParticleTypeDatabase::Proton.Mass());
            const LorentzVec missing = beam_target - photon_sum;
            const double missing_mass = missing.M();

            h_missingmass.Fill(missing_mass);
            const interval<double> MM_cut(800, 1050);
            if(!MM_cut.Contains(missing_mass))
                continue;

            steps->Fill("Fitted taggerhits",1.0);

            // more sophisticated fitter
            fitter.SetEgammaBeam(taggerhit.PhotonEnergy);
            fitter.SetProton(proton);
            fitter.SetPhotons(photons);
            auto fit_result = fitter.DoFit();

            if(fit_result.Status == APLCON::Result_Status_t::Success) {
                steps->Fill("Fit OK",1.0);
                h_fitprobability.Fill(fit_result.Probability);
                utils::ParticleTools::FillIMCombinations([this] (double x) {IM_2g_byFit.Fill(x);},  2, photons);
                utils::ParticleTools::FillIMCombinations([this] (double x) {IM_2g_fitted.Fill(x);},  2, fitter.GetFittedPhotons());
                FillPulls(fit_result.Variables);
                tree->Fill();
            }
            else {
                steps->Fill("Fit failed",1.0);
            }

        }
    }
}

void KinFitPi0::MultiPi0::FillPulls(const fit_variables_t& vars)
{
    for(auto& var : vars) {
        const auto& varname = var.first;
        const auto& pull = var.second.Pull;
        // get the iterator to where the key *should* be if it existed:
        auto hint = h_pulls.lower_bound(varname);
        if (hint == h_pulls.end() || h_pulls.key_comp()(varname, hint->first)) {
            // insert at the correct location
            hint = h_pulls.emplace_hint(hint, varname, PromptRandom::Hist1(promptrandom));
            hint->second.MakeHistograms(HistFac,"pull_"+varname,varname+": Pull",BinSettings(100,-3,3),"Pull","");
        }
        hint->second.Fill(pull);
    }
}

void KinFitPi0::MultiPi0::ShowResult()
{
    canvas(std_ext::formatter() << "KinFitPi0: " << multiplicity << "Pi0")
            << steps
            << Proton_Coplanarity
            << h_taggtime
            << h_missingmass.subtracted
            << h_fitprobability.subtracted
            << IM_2g_byFit.subtracted
            << IM_2g_fitted.subtracted
            << endc;

    canvas c_pulls(std_ext::formatter() << "KinFitPi0: " << multiplicity << "Pulls");

    for(auto& h_pull : h_pulls)
        c_pulls << h_pull.second.subtracted;

    c_pulls << endc;
}

KinFitPi0::MultiPi0::IM_perms_t KinFitPi0::MultiPi0::BuildIMPerms(unsigned multiplicity)
{
    ParticleTypeTree t;
    if(multiplicity==1)
        t = ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0_2g);
    else if(multiplicity==2)
        t = ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::TwoPi0_4g);
    else if(multiplicity==3)
        t = ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::ThreePi0_6g);
    if(!t)
        throw runtime_error("Multiplicity not supported");
    auto t_tmp = t->DeepCopy<int>([] (const ParticleTypeTree&) {return 0;});
    for(const auto& daughter : t_tmp->Daughters()) {
        if(daughter->IsLeaf()) {
            daughter->Unlink();
            break;
        }
    }
    t_tmp->Sort();
    IM_perms_t perms;
    using t_tmp_t = decltype(t_tmp);
    std::vector<t_tmp_t::element_type::node_t> leaves;
    t_tmp->GetUniquePermutations(leaves, perms);
    return perms;
}



AUTO_REGISTER_PHYSICS(KinFitPi0)
