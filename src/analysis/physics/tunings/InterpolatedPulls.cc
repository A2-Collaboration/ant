#include "InterpolatedPulls.h"

#include "base/std_ext/misc.h"
#include "base/Logger.h"
#include "expconfig/ExpConfig.h"
#include "utils/ParticleTools.h"
#include "analysis/utils/uncertainties/FitterSergey.h"
#include "utils/ProtonPhotonCombs.h"

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace std;

InterpolatedPulls::InterpolatedPulls(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    promptrandom(ExpConfig::Setup::Get()),
    fit_model_data(utils::UncertaintyModels::Interpolated::makeAndLoad(
                   utils::UncertaintyModels::Interpolated::Type_t::Data,
                   // use Sergey as starting point
                   make_shared<utils::UncertaintyModels::FitterSergey>()
                  )
          ),
    fit_model_mc(utils::UncertaintyModels::Interpolated::makeAndLoad(
                 utils::UncertaintyModels::Interpolated::Type_t::MC,
                 // use Sergey as starting point
                 make_shared<utils::UncertaintyModels::FitterSergey>()
                )
          ),
    fitter(nullptr, // fit model will be set event-by-event
           opts->Get<bool>("FitZVertex", true) // enable Z vertex by default
           ),
    pullswriter(HistFac)
{

    if(fitter.IsZVertexFitEnabled()) {
        fitter.SetZVertexSigma(0); // use unmeasured z vertex
        LOG(INFO) << "Running with unmeasured Z Vertex";
    } else {
        LOG(INFO) << "Running with fixed Z Vertex";
    }

    steps = HistFac.makeTH1D("Steps","","#",BinSettings(10),"steps");

    BinSettings bins_MM(200,640,1240);
    BinSettings bins_E(300,0,1000);
    BinSettings bins_vetoE(100,0,10);
    BinSettings bins_ToF(100,-10,15);
    BinSettings bins_PSA_phi(100,25,65);
    BinSettings bins_IM_gg(300,0,750);

    h_missingmass_cut = HistFac.makeTH1D("MissingMass","MM / MeV","",
                                         bins_MM,"h_missingmass_cut");
    h_missingmass_best = HistFac.makeTH1D("MissingMass","MM / MeV","",
                                          bins_MM,"h_missingmass_best");

    h_IM_gg_gg = HistFac.makeTH2D("IM 2#gamma 2#gamma","IM / MeV","IM / MeV",bins_IM_gg, bins_IM_gg,"h_IM_gg_gg");
    h_IM_gg_gg_cut = HistFac.makeTH2D("IM 2#gamma 2#gamma","IM / MeV","IM / MeV",bins_IM_gg, bins_IM_gg,"h_IM_gg_gg_cut");

    h_zvertex = HistFac.makeTH1D("z Vertex","z / cm","",
                                 BinSettings(50,-15,15),
                                 "h_zvertex");

    h_proton_E_theta = HistFac.makeTH2D("Proton","E / MeV","#theta / #circ",
                                        bins_E,BinSettings(200,0,100),"h_proton_E_theta");

    h_ToF_E_photon_taps = HistFac.makeTH2D("Photon TAPS ToF","#Delta t / ns","E / MeV",
                                           bins_ToF,bins_E,"h_ToF_E_photon_taps");
    h_ToF_E_proton_taps = HistFac.makeTH2D("Proton TAPS ToF","#Delta t / ns","E / MeV",
                                           bins_ToF,bins_E,"h_ToF_E_proton_taps");

    h_PSA_photon_taps = HistFac.makeTH2D("Photon TAPS PSA","#phi / #circ","R / MeV",
                                         bins_PSA_phi,bins_E,"h_PSA_photon_taps");
    h_PSA_proton_taps = HistFac.makeTH2D("Proton TAPS PSA","#phi / #circ","R / MeV",
                                         bins_PSA_phi,bins_E,"h_PSA_proton_taps");

    h_E_vetoE_photon_cb = HistFac.makeTH2D("Photon CB dE/E","E / MeV","#Delta E / MeV",
                                           bins_E,bins_vetoE,"h_E_vetoE_photon_cb");
    h_E_vetoE_photon_taps = HistFac.makeTH2D("Photon TAPS dE/E","E / MeV","#Delta E / MeV",
                                             bins_E,bins_vetoE,"h_E_vetoE_photon_taps");
    h_E_vetoE_proton_cb = HistFac.makeTH2D("Proton CB dE/E","E / MeV","#Delta E / MeV",
                                           bins_E,bins_vetoE,"h_E_vetoE_proton_cb");
    h_E_vetoE_proton_taps = HistFac.makeTH2D("Proton TAPS dE/E","E / MeV","#Delta E / MeV",
                                             bins_E,bins_vetoE,"h_E_vetoE_proton_taps");

    t.CreateBranches(HistFac.makeTTree("t"));
}

void InterpolatedPulls::ProcessEvent(const TEvent& event, manager_t&)
{
    triggersimu.ProcessEvent(event);

    const TEventData& data = event.Reconstructed();

    steps->Fill("Seen",1);

    // cut on energy sum and number of candidates
    if(!triggersimu.HasTriggered())
        return;
    steps->Fill("Triggered",1);

    const auto& cands = data.Candidates;
    if(cands.size() != 5)
        return;
    steps->Fill("nCands==5",1);

    // choose uncertainty depending on Data/MC input
    const bool is_MC = data.ID.isSet(TID::Flags_t::MC);
    fitter.SetUncertaintyModel(is_MC ? fit_model_mc : fit_model_data);

    utils::ProtonPhotonCombs proton_photons(data.Candidates);

    for(const TTaggerHit& taggerhit : data.TaggerHits) {

        steps->Fill("Seen taggerhits",1.0);

        promptrandom.SetTaggerTime(triggersimu.GetCorrectedTaggerTime(taggerhit));
        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        const auto& TaggW = promptrandom.FillWeight();

        auto filtered_combs = proton_photons()
                              .Observe([this] (const string& cut) { steps->Fill(cut.c_str(), 1.0); }, "F ")
                              .FilterMM(taggerhit, ParticleTypeDatabase::Proton.GetWindow(300).Round());

        if(filtered_combs.empty()) {
            steps->Fill("No combs left",1.0);
            continue;
        }

        double best_prob = std_ext::NaN;
        std::vector<utils::Fitter::FitParticle> best_fitParticles;
        double best_zvertex = std_ext::NaN;
        t.IM_pi0pi0_radius = std_ext::NaN;

        // use any candidate as proton, and do the analysis (ignore ParticleID stuff)

        for(const utils::ProtonPhotonCombs::comb_t& comb : filtered_combs) {

            h_missingmass_cut->Fill(comb.MissingMass, TaggW);

            auto IM_pi0pi0_radius = std_ext::NaN;
            vec2 IM_pi0pi0_best;

            const vector<vector<unsigned>> goldhaber_comb{{0,1,2,3},{0,2,1,3},{0,3,1,2}};
            for(auto& i : goldhaber_comb) {
                const auto& p = comb.Photons;

                const vec2 IM_pi0pi0_true(ParticleTypeDatabase::Pi0.Mass(), ParticleTypeDatabase::Pi0.Mass());
                const vec2 IM_pi0pi0(
                            (*p[i[0]] + *p[i[1]]).M(),
                            (*p[i[2]] + *p[i[3]]).M());

                h_IM_gg_gg->Fill(IM_pi0pi0.x, IM_pi0pi0.y, TaggW);

               if(std_ext::copy_if_better(IM_pi0pi0_radius, (IM_pi0pi0 - IM_pi0pi0_true).R(), std::less<double>())) {
                   IM_pi0pi0_best = IM_pi0pi0;
               }
            }

            // this cut is crucial to include include only pi0pi0 events
            if(IM_pi0pi0_radius > 15)
                continue;
            steps->Fill("#pi^{0}#pi^{0}",1.0);

            for(auto& i : goldhaber_comb) {
                const auto& p = comb.Photons;
                const auto& IM1 = (*p[i[0]] + *p[i[1]]).M();
                const auto& IM2 = (*p[i[2]] + *p[i[3]]).M();

                h_IM_gg_gg_cut->Fill(IM1, IM2, TaggW);
            }

            // do the fitting

            const auto& fit_result = fitter.DoFit(taggerhit.PhotonEnergy, comb.Proton, comb.Photons);

            if(fit_result.Status != APLCON::Result_Status_t::Success)
                continue;

            steps->Fill("KinFit OK",1);

            if(!std_ext::copy_if_greater(best_prob, fit_result.Probability))
                continue;

            best_fitParticles = fitter.GetFitParticles();
            best_zvertex = fitter.GetFittedZVertex();
            t.IM_pi0pi0_radius = IM_pi0pi0_radius;
            t.IM_gg_gg().at(0) = IM_pi0pi0_best.x;
            t.IM_gg_gg().at(1) = IM_pi0pi0_best.y;
        }

        if(!isfinite(best_prob))
            continue;

        steps->Fill("Fill",1);

        h_zvertex->Fill(best_zvertex, TaggW);
        pullswriter.Fill(best_fitParticles, TaggW, best_prob, best_zvertex, t);

        // fill the many check hists
        LorentzVec best_photon_sum({0,0,0},0);
        for(const utils::Fitter::FitParticle& fitparticle : best_fitParticles) {
            const TCandidatePtr& cand = fitparticle.Particle->Candidate;
            const TParticlePtr& fitted = fitparticle.AsFitted();

            const vec2 PSA{fitted->Ek(), cand->FindCaloCluster()->ShortEnergy};

            if(fitparticle.Particle->Type() == ParticleTypeDatabase::Photon) {
                // photon
                best_photon_sum += *fitparticle.Particle;
                if(cand->Detector & Detector_t::Type_t::CB) {
                    // CB
                    h_E_vetoE_photon_cb->Fill(fitted->Ek(), cand->VetoEnergy);
                }
                else if(cand->Detector & Detector_t::Type_t::TAPS) {
                    // TAPS
                    h_ToF_E_photon_taps->Fill(cand->Time-triggersimu.GetRefTiming(), fitted->Ek(), TaggW);
                    if(fitted->Ek()<400)
                        h_PSA_photon_taps->Fill(std_ext::radian_to_degree(PSA.Phi()), PSA.R());
                    h_E_vetoE_photon_taps->Fill(fitted->Ek(), cand->VetoEnergy);

                }
            }
            else  if(fitparticle.Particle->Type() == ParticleTypeDatabase::Proton) {
                // proton
                h_proton_E_theta->Fill(fitted->Ek(), std_ext::radian_to_degree(fitted->Theta()), TaggW);
                if(cand->Detector & Detector_t::Type_t::CB) {
                    // CB
                    h_E_vetoE_proton_cb->Fill(fitted->Ek(), cand->VetoEnergy);
                }
                else if(cand->Detector & Detector_t::Type_t::TAPS) {
                    // TAPS
                    h_ToF_E_proton_taps->Fill(cand->Time-triggersimu.GetRefTiming(), fitted->Ek(), TaggW);
                    if(fitted->Ek()<400)
                        h_PSA_proton_taps->Fill(std_ext::radian_to_degree(PSA.Phi()), PSA.R());
                    h_E_vetoE_proton_taps->Fill(fitted->Ek(), cand->VetoEnergy);
                }
            }
        }

        if(best_prob>0.01) {
            const LorentzVec& best_missing = taggerhit.GetPhotonBeam() + LorentzVec({0, 0, 0}, ParticleTypeDatabase::Proton.Mass()) - best_photon_sum;
            h_missingmass_best->Fill(best_missing.M(), TaggW);
        }

        // check some MCTrue stuff if available
        auto& particletree = event.MCTrue().ParticleTree;
        if(particletree) {
            steps->Fill("MCTrue",1);
            // 1=Signal, 2=Reference, 9=MissedBkg, >=10 found in ptreeBackgrounds
            if(particletree->IsEqual(ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::TwoPi0_4g),
                                     utils::ParticleTools::MatchByParticleName)) {
                steps->Fill("MC #pi^{0}#pi^{0}",1);
            }
            else {
                steps->Fill("MC Bkg",1);
            }
        }
    }
}

void InterpolatedPulls::ShowResult()
{
    canvas("Overview") << steps << h_missingmass_best
                       << h_zvertex
                       << drawoption("colz")
                       << h_IM_gg_gg_cut
                       << endc;
}

void InterpolatedPulls::Finish()
{
    LOG(INFO) << "Fit Model Statistics Data:\n" << *fit_model_data;
    LOG(INFO) << "Fit Model Statistics MC:\n" << *fit_model_mc;
}

AUTO_REGISTER_PHYSICS(InterpolatedPulls)
