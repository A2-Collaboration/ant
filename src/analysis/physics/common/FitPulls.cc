#include "FitPulls.h"
#include "base/std_ext/vector.h"
#include "TH1D.h"

#include "utils/particle_tools.h"
#include "plot/root_draw.h"

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace std;

const vector<ParticleTypeTreeDatabase::Channel> FitPulls::channels = {
    ParticleTypeTreeDatabase::Channel::Pi0_2g,
    ParticleTypeTreeDatabase::Channel::TwoPi0_4g,
    ParticleTypeTreeDatabase::Channel::ThreePi0_6g,
    ParticleTypeTreeDatabase::Channel::Omega_gPi0_3g,
};

FitPulls::FitPulls(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    opt_save_after_cut(opts->Get<bool>("SaveAfterCut", false)),
    opt_save_only(opts->Get<bool>("SaveOnly", false))
{
    promptrandom.AddPromptRange({-2.5,1.5}); // slight offset due to CBAvgTime reference
    promptrandom.AddRandomRange({-50,-10});  // just ensure to be way off prompt peak
    promptrandom.AddRandomRange({  10,50});

    auto uncertainty_model = make_shared<utils::UncertaintyModels::Optimized_Oli1>();

    const auto& histFac_name = uncertainty_model->to_string_simple();
    HistogramFactory histFac(histFac_name, HistFac);

    h_protoncopl = histFac.makeTH1D("Coplanarity","#delta#phi / degree","",BinSettings(100,-180,180),"h_protoncopl");
    h_taggtime = histFac.makeTH1D("Tagged Time","t / ns", "", BinSettings(300,-60,60), "h_taggtime");

    auto make_hstack = [histFac_name, histFac] (const std::string& name) {
        return histFac.make<ant::hstack>(name, "hstack: "+histFac_name+": "+name);
    };

    hstacks.emplace_back(make_hstack("h_probability"));
    hstacks.emplace_back(make_hstack("p_cb_g_E"));
    hstacks.emplace_back(make_hstack("p_cb_g_Phi"));
    hstacks.emplace_back(make_hstack("p_cb_g_Theta"));
    hstacks.emplace_back(make_hstack("p_cb_p_E"));
    hstacks.emplace_back(make_hstack("p_cb_p_Theta"));
    hstacks.emplace_back(make_hstack("p_cb_p_Phi"));
    hstacks.emplace_back(make_hstack("p_taps_g_E"));
    hstacks.emplace_back(make_hstack("p_taps_g_Theta"));
    hstacks.emplace_back(make_hstack("p_taps_g_Phi"));
    hstacks.emplace_back(make_hstack("p_taps_p_E"));
    hstacks.emplace_back(make_hstack("p_taps_p_Theta"));
    hstacks.emplace_back(make_hstack("p_taps_p_Phi"));


    // create fitter for each channel
    unsigned n_ch = 0;
    for(auto& channel : channels) {
        ChannelItem_t item(histFac, channel, uncertainty_model);
        ChannelHists_t& hists = *item.Hists;

        std::vector<TH1*> histptrs{
                    hists.h_probability,
                    hists.p_cb_g_E,
                    hists.p_cb_g_Theta,
                    hists.p_cb_g_Phi,
                    hists.p_cb_p_E,
                    hists.p_cb_p_Theta,
                    hists.p_cb_p_Phi,
                    hists.p_taps_g_E,
                    hists.p_taps_g_Theta,
                    hists.p_taps_g_Phi,
                    hists.p_taps_p_E,
                    hists.p_taps_p_Theta,
                    hists.p_taps_p_Phi
        };

        for(unsigned i=0;i<histptrs.size();i++) {
            histptrs[i]->SetLineColor(ColorPalette::Colors[n_ch]);
            *hstacks.at(i) << histptrs[i];
        }

        auto& items = treefitters[item.Multiplicity];
        items.emplace_back(move(item));
        n_ch++;
    }


}

FitPulls::ChannelItem_t::ChannelItem_t(const HistogramFactory& parent,
                                     ParticleTypeTreeDatabase::Channel channel,
                                     utils::Fitter::UncertaintyModelPtr uncertainty_model
                                     )
{
    auto ptree = ParticleTypeTreeDatabase::Get(channel);
    auto count_leaves = [] (const ParticleTypeTree& tree) {
        unsigned leaves = 0;
        tree->Map_nodes(
                    [&leaves] (const ParticleTypeTree& n) { if(n->IsLeaf()) leaves++; }
        );
        return leaves;
    };
    Multiplicity = count_leaves(ptree);

    const auto& ptree_name = utils::ParticleTools::GetDecayString(ptree);

    Fitter = std_ext::make_unique<utils::TreeFitter>("treefitter_" + ptree_name,
                                                     ptree,
                                                     Multiplicity-1, // assume one particle to be proton
                                                     uncertainty_model);

    Hists = std_ext::make_unique<ChannelHists_t>(parent, ptree_name);
}

void FitPulls::findProton(const TCandidateList& cands,
                          const TTaggerHit& taggerhit,
                          TParticlePtr& proton,
                          TParticleList& photons,
                          LorentzVec& photon_sum)
{
    // find the candidate which gives the closest missing mass to proton
    auto min_mm = std::numeric_limits<double>::infinity();
    for(auto& cand_proton : cands.get_iter()) {
        LorentzVec photon_sum_tmp(0,0,0,0);
        for(auto& cand : cands.get_iter()) {
            if(cand == cand_proton)
                continue;
            photon_sum_tmp += TParticle(ParticleTypeDatabase::Photon, cand.get_ptr());
        }
        const LorentzVec beam_target = taggerhit.GetPhotonBeam() + LorentzVec(0, 0, 0, ParticleTypeDatabase::Proton.Mass());
        const LorentzVec missing = beam_target - photon_sum_tmp;
        const double mm = abs(missing.M() - ParticleTypeDatabase::Proton.Mass());
        if(mm < min_mm) {
            min_mm = mm;
            proton = make_shared<TParticle>(ParticleTypeDatabase::Proton, cand_proton.get_ptr());
            photon_sum = photon_sum_tmp;
        }
    }

    // fill the photons
    for(auto& cand_photon : cands.get_iter()) {
        auto cand_ptr = cand_photon.get_ptr();
        if(cand_ptr == proton->Candidate)
            continue;
        photons.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::Photon, cand_ptr));
    }
}

void FitPulls::ProcessEvent(const TEvent& event, manager_t& manager)
{
    const auto& cands = event.Reconstructed().Candidates;

    auto fitters = treefitters.find(unsigned(cands.size()));
    if(fitters == treefitters.end())
        return;

    // check that all candidates are clean
    for(const TCandidate& cand : cands) {
        if(auto calocluster = cand.FindCaloCluster()) {
            if(calocluster->HasFlag(TCluster::Flags_t::TouchesHole))
                return;
        }
        else
            return;
    }

    if(opt_save_after_cut) {
        manager.SaveEvent();
        if(opt_save_only)
            return;
    }

    for(const TTaggerHit& taggerhit : event.Reconstructed().TaggerHits) {
        const auto taggtime = taggerhit.Time - event.Reconstructed().Trigger.CBTiming;
        promptrandom.SetTaggerHit(taggtime);
        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        h_taggtime->Fill(taggtime);

        // select the proton as the particle with best missing mass
        TParticlePtr  proton;
        TParticleList photons;
        LorentzVec photon_sum;
        findProton(cands, taggerhit, proton, photons, photon_sum);

        // check coplanarity
        const double d_phi = std_ext::radian_to_degree(vec2::Phi_mpi_pi(proton->Phi()-photon_sum.Phi() - M_PI ));
        h_protoncopl->Fill(d_phi, promptrandom.FillWeight());
        const interval<double> copl_cut{-20, 20};
        if(!copl_cut.Contains(d_phi))
            continue;

        // do all the fits
        for(const ChannelItem_t& item : fitters->second) {
            utils::TreeFitter& fitter = *item.Fitter;
            fitter.SetEgammaBeam(taggerhit.PhotonEnergy);
            fitter.SetPhotons(photons);
            fitter.SetProton(proton);

            APLCON::Result_t result;
            double max_prob = 0;
            std::vector<utils::Fitter::FitParticle> fitparticles;
            while(fitter.NextFit(result)) {
                if(result.Status != APLCON::Result_Status_t::Success)
                    continue;
                if(result.Probability<max_prob)
                    continue;
                max_prob = result.Probability;
                fitparticles = fitter.GetFitParticles();
            }

            const ChannelHists_t& h = *item.Hists;

            h.h_probability->Fill(max_prob, promptrandom.FillWeight());

            for(const auto& fitparticle : fitparticles) {
                const auto& p = fitparticle.Particle;
                // select the right set of histograms
                auto& h_E = p->Type() == ParticleTypeDatabase::Photon ?
                                (p->Candidate->Detector & Detector_t::Type_t::CB ? h.p_cb_g_E : h.p_taps_g_E)
                              :
                                (p->Candidate->Detector & Detector_t::Type_t::CB ? h.p_cb_p_E : h.p_taps_p_E);
                h_E->Fill(fitparticle.Ek.Pull, promptrandom.FillWeight());

                auto& h_Theta = p->Type() == ParticleTypeDatabase::Photon ?
                                    (p->Candidate->Detector & Detector_t::Type_t::CB ? h.p_cb_g_Theta : h.p_taps_g_Theta)
                                  :
                                    (p->Candidate->Detector & Detector_t::Type_t::CB ? h.p_cb_p_Theta : h.p_taps_p_Theta);
                h_Theta->Fill(fitparticle.Theta.Pull, promptrandom.FillWeight());

                auto& h_Phi = p->Type() == ParticleTypeDatabase::Photon ?
                                  (p->Candidate->Detector & Detector_t::Type_t::CB ? h.p_cb_g_Phi : h.p_taps_g_Phi)
                                :
                                  (p->Candidate->Detector & Detector_t::Type_t::CB ? h.p_cb_p_Phi : h.p_taps_p_Phi);
                h_Phi->Fill(fitparticle.Phi.Pull, promptrandom.FillWeight());
            }
        }
    }
}

void FitPulls::ShowResult()
{

}





FitPulls::ChannelHists_t::ChannelHists_t(const HistogramFactory& h, const string& name)
{
    HistogramFactory histFac(name, h, name);

    h_probability = histFac.makeTH1D("Probability","p","",BinSettings(100,0,1),"h_probability");

    BinSettings bins_pulls(30,-3,3);
    p_cb_g_E     = histFac.makeTH1D("p_cb_g_E",    "","",bins_pulls,"p_cb_g_E");
    p_cb_g_Theta = histFac.makeTH1D("p_cb_g_Theta","","",bins_pulls,"p_cb_g_Theta");
    p_cb_g_Phi   = histFac.makeTH1D("p_cb_g_Phi",  "","",bins_pulls,"p_cb_g_Phi");
    p_cb_p_E     = histFac.makeTH1D("p_cb_p_E",    "","",bins_pulls,"p_cb_p_E");
    p_cb_p_Theta = histFac.makeTH1D("p_cb_p_Theta","","",bins_pulls,"p_cb_p_Theta");
    p_cb_p_Phi   = histFac.makeTH1D("p_cb_p_Phi",  "","",bins_pulls,"p_cb_p_Phi");
    p_taps_g_E     = histFac.makeTH1D("p_taps_g_E",    "","",bins_pulls,"p_taps_g_E");
    p_taps_g_Theta = histFac.makeTH1D("p_taps_g_Theta","","",bins_pulls,"p_taps_g_Theta");
    p_taps_g_Phi   = histFac.makeTH1D("p_taps_g_Phi",  "","",bins_pulls,"p_taps_g_Phi");
    p_taps_p_E     = histFac.makeTH1D("p_taps_p_E",    "","",bins_pulls,"p_taps_p_E");
    p_taps_p_Theta = histFac.makeTH1D("p_taps_p_Theta","","",bins_pulls,"p_taps_p_Theta");
    p_taps_p_Phi   = histFac.makeTH1D("p_taps_p_Phi",  "","",bins_pulls,"p_taps_p_Phi");
}





AUTO_REGISTER_PHYSICS(FitPulls)
