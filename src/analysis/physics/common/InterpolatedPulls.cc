#include "InterpolatedPulls.h"

#include "plot/root_draw.h"
#include "base/std_ext/misc.h"

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace std;

InterpolatedPulls::InterpolatedPulls(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    model(utils::UncertaintyModels::Interpolated::makeAndLoad(
              // use OptimizedOli1 as default
              make_shared<utils::UncertaintyModels::Optimized_Oli1>()
              )
          ),
    fitter("KinFit", 4, model, true),
    pullswriter(HistFac)
{
    fitter.SetZVertexSigma(0);


    promptrandom.AddPromptRange({-2.5,2.5});
    promptrandom.AddRandomRange({-50,-5});
    promptrandom.AddRandomRange({  5,50});

    steps = HistFac.makeTH1D("Steps","","#",BinSettings(15),"steps");

}

void InterpolatedPulls::ProcessEvent(const TEvent& event, manager_t&)
{
    const TEventData& data = event.Reconstructed();

    steps->Fill("Seen",1);

    // cut on energy sum and number of candidates
    if(data.Trigger.CBEnergySum <= 550)
        return;
    steps->Fill("CBESum>550MeV",1);

    const auto& cands = data.Candidates;
    if(cands.size() != 5)
        return;
    steps->Fill("nCands==5",1);



    for(const TTaggerHit& taggerhit : data.TaggerHits) {


        steps->Fill("Seen taggerhits",1.0);

        promptrandom.SetTaggerHit(taggerhit.Time);
        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;

//        TParticlePtr  selected_proton;
//        TParticleList selected_photons;
        std::vector<utils::Fitter::FitParticle> best_fitParticles;
        double best_prob = std_ext::NaN;

        // use any candidate as proton, and do the analysis (ignore ParticleID stuff)
        for(auto i_proton : cands.get_iter()) {

            steps->Fill("Seen protons",1.0);

            const auto proton = std::make_shared<TParticle>(ParticleTypeDatabase::Proton, i_proton);
            std::vector<TParticlePtr> photons;
            for(auto i_photon : cands.get_iter()) {
                if(i_photon == i_proton)
                    continue;
                photons.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::Photon, i_photon));
            }


            LorentzVec photon_sum(0,0,0,0);
            for(const auto& p : photons) {
                photon_sum += *p;
            }

            // proton coplanarity
            const double d_phi = std_ext::radian_to_degree(vec2::Phi_mpi_pi(proton->Phi()-photon_sum.Phi() - M_PI ));
            if(d_phi<-20 || d_phi>20)
                continue;
            steps->Fill("Copl p in [-20;20]",1);

            // missing mass
            const LorentzVec beam_target = taggerhit.GetPhotonBeam() + LorentzVec(0, 0, 0, ParticleTypeDatabase::Proton.Mass());
            const LorentzVec missing = beam_target - photon_sum;
            const double missing_mass = missing.M();

            if(missing_mass<840 || missing_mass>1040)
                continue;
            steps->Fill("MM in [840;1040]",1);

            auto angle_p_calcp = std_ext::radian_to_degree(missing.Angle(proton->p));
            if(angle_p_calcp > 15.0)
                continue;
            steps->Fill("p angle < 15.0#circ",1);

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



            fitter.SetEgammaBeam(taggerhit.PhotonEnergy);
            fitter.SetProton(proton);
            fitter.SetPhotons(photons);
            const auto& fit_result = fitter.DoFit();


            if(fit_result.Status != APLCON::Result_Status_t::Success)
                continue;

            steps->Fill("KinFit OK",1);

            if(!std_ext::copy_if_greater(best_prob, fit_result.Probability))
                continue;

            best_fitParticles = fitter.GetFitParticles();
        }

        if(!isfinite(best_prob))
            continue;

        steps->Fill("Fill",1);

        pullswriter.Fill(best_fitParticles, promptrandom.FillWeight(), best_prob);
    }

}

void InterpolatedPulls::ShowResult()
{
    canvas("Overview") << steps << endc;
}

void InterpolatedPulls::Finish()
{

}

AUTO_REGISTER_PHYSICS(InterpolatedPulls)
