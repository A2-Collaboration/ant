#include "InterpolatedPulls.h"

#include "base/std_ext/misc.h"
#include "base/Logger.h"
#include "expconfig/ExpConfig.h"
#include "utils/ParticleTools.h"
#include "analysis/utils/uncertainties/FitterSergey.h"
#include "analysis/utils/uncertainties/Optimized.h"

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace std;

InterpolatedPulls::InterpolatedPulls(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    promptrandom(ExpConfig::Setup::Get()),
    fit_model(utils::UncertaintyModels::Interpolated::makeAndLoad(
                  // use Sergey as starting point
                  make_shared<utils::UncertaintyModels::FitterSergey>()
                  )
          ),
    fitter(fit_model,
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

    h_missingmass = HistFac.makeTH1D("MissingMass","MM / MeV","",
                                     bins_MM,"h_missingmass");
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

    for(const TTaggerHit& taggerhit : data.TaggerHits) {

        steps->Fill("Seen taggerhits",1.0);

        promptrandom.SetTaggerTime(triggersimu.GetCorrectedTaggerTime(taggerhit));
        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        const auto& TaggW = promptrandom.FillWeight();

        double best_prob = std_ext::NaN;
        std::vector<utils::Fitter::FitParticle> best_fitParticles;
        double best_zvertex = std_ext::NaN;

        // use any candidate as proton, and do the analysis (ignore ParticleID stuff)
        for(auto i_proton : cands.get_iter()) {

            steps->Fill("Seen protons",1.0);

            TParticlePtr proton = std::make_shared<TParticle>(ParticleTypeDatabase::Proton, i_proton);
            std::vector<TParticlePtr> photons;
            for(auto i_photon : cands.get_iter()) {
                if(i_photon == i_proton)
                    continue;
                photons.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::Photon, i_photon));
            }


            LorentzVec photon_sum({0,0,0},0);
            for(const auto& p : photons) {
                photon_sum += *p;
            }

            // proton coplanarity
            const double d_phi = std_ext::radian_to_degree(vec2::Phi_mpi_pi(proton->Phi()-photon_sum.Phi() - M_PI ));
            if(d_phi<-20 || d_phi>20)
                continue;
            steps->Fill("Copl p in [-20;20]",1);

            // missing mass
            const LorentzVec& beam_target = taggerhit.GetPhotonBeam() + LorentzVec({0, 0, 0}, ParticleTypeDatabase::Proton.Mass());
            const LorentzVec& missing = beam_target - photon_sum;
            const double missing_mass = missing.M();
            h_missingmass->Fill(missing_mass, TaggW);

            if(missing_mass<840 || missing_mass>1040)
                continue;
            steps->Fill("MM in [840;1040]",1);

            auto angle_p_calcp = std_ext::radian_to_degree(missing.Angle(proton->p));
            if(angle_p_calcp > 15.0)
                continue;
            steps->Fill("p angle < 15.0#circ",1);

            h_missingmass_cut->Fill(missing_mass, TaggW);

            // check gammas
            bool is_Pi0Pi0 = false;
            bool is_Pi0Eta = false;
            const auto& Pi0_cut = ParticleTypeDatabase::Pi0.GetWindow(30);
            const auto& Eta_cut = ParticleTypeDatabase::Eta.GetWindow(60);

            const vector<vector<unsigned>> goldhaber_comb{{0,1,2,3},{0,2,1,3},{0,3,1,2}};
            for(auto& i : goldhaber_comb) {
                const auto& p = photons;
                const auto& IM1 = (*p[i[0]] + *p[i[1]]).M();
                const auto& IM2 = (*p[i[2]] + *p[i[3]]).M();

                h_IM_gg_gg->Fill(IM1, IM2, TaggW);

                if(Pi0_cut.Contains(IM1) && Pi0_cut.Contains(IM2))
                    is_Pi0Pi0 = true;
                if(   (Pi0_cut.Contains(IM1) && Eta_cut.Contains(IM2))
                   || (Pi0_cut.Contains(IM2) && Eta_cut.Contains(IM1))
                  )
                    is_Pi0Eta = true;
            }

            if(!is_Pi0Pi0 && !is_Pi0Eta)
                continue;
            steps->Fill("#pi^{0}#pi^{0}",is_Pi0Pi0);
            steps->Fill("#pi^{0}#eta",is_Pi0Eta);

            for(auto& i : goldhaber_comb) {
                const auto& p = photons;
                const auto& IM1 = (*p[i[0]] + *p[i[1]]).M();
                const auto& IM2 = (*p[i[2]] + *p[i[3]]).M();

                h_IM_gg_gg_cut->Fill(IM1, IM2, TaggW);
            }

            // do the fitting

            const auto& fit_result = fitter.DoFit(taggerhit.PhotonEnergy, proton, photons);

            if(fit_result.Status != APLCON::Result_Status_t::Success)
                continue;

            steps->Fill("KinFit OK",1);

            if(!std_ext::copy_if_greater(best_prob, fit_result.Probability))
                continue;

            best_fitParticles = fitter.GetFitParticles();
            best_zvertex = fitter.GetFittedZVertex();
        }

        if(!isfinite(best_prob))
            continue;

        steps->Fill("Fill",1);

        h_zvertex->Fill(best_zvertex, TaggW);

        pullswriter.Fill(best_fitParticles, TaggW, best_prob, best_zvertex);

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
            else if(particletree->IsEqual(ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0Eta_4g),
                                          utils::ParticleTools::MatchByParticleName)) {
                steps->Fill("MC #pi^{0}#eta",1);
            }
            else {
                steps->Fill("MC Bkg",1);
            }
        }
    }

}

void InterpolatedPulls::ShowResult()
{
    canvas("Overview") << steps << h_missingmass_best << endc;
}

void InterpolatedPulls::Finish()
{
    if(fit_model) {
        LOG(INFO) << "Fit Model Statistics:\n" << *fit_model;
    }
}

AUTO_REGISTER_PHYSICS(InterpolatedPulls)
