#include "JustParticles.h"
#include "base/std_ext/math.h"
#include "expconfig/ExpConfig.h"
#include "analysis/utils/uncertainties/FitterSergey.h"

#include "TTree.h"

#include <cassert>

using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace std;


JustParticles::JustParticles(const string& name, OptionsPtr opts):
    Physics(name, opts),
    multiplicities(opts->Get<decltype(multiplicities)>("PhotonMulti",{{2,2},{4,4},{6,6}})),
    save_events(opts->Get<bool>("SaveEvents",false))
{
    auto enclosing = multiplicities.EnclosingInterval();
    if(!enclosing.IsSane() || enclosing.Start()<1)
        throw runtime_error("Given photon multiplicities not sane");


    steps = HistFac.makeTH1D("Steps","","",BinSettings(10),"steps");

    promptrandom.AddPromptRange({-5, 5}); // slight offset due to CBAvgTime reference
    promptrandom.AddRandomRange({-50,-10});  // just ensure to be way off prompt peak
    promptrandom.AddRandomRange({  10, 50});

    tree = HistFac.makeTTree("tree");
    t.CreateBranches(tree);

    // prepare fitters for all multiplicities
    fitters.resize(enclosing.Stop());

    auto uncertainty = make_shared<utils::UncertaintyModels::FitterSergey>();

    for(unsigned mult=enclosing.Start();mult<=enclosing.Stop();mult++) {
        if(!multiplicities.Contains(mult))
            continue;
        auto fitter = std_ext::make_unique<utils::KinFitter>(uncertainty);

        fitters[mult-1] = move(fitter);
    }

    // needed for calculating ToF
    taps_detector = ExpConfig::Setup::GetDetector<expconfig::detector::TAPS>();
}

void JustParticles::ProcessEvent(const TEvent& event, manager_t& manager)
{
    triggersimu.ProcessEvent(event);

    steps->Fill("Seen",1.0);

    const auto& cands = event.Reconstructed().Candidates;
    TCandidatePtrList cands_taps;
    TCandidatePtrList cands_cb;

    t.b_CBSumVetoE = 0;
    for(const auto& p : cands.get_iter()) {
        if(p->Detector & Detector_t::Any_t::TAPS_Apparatus) {
            cands_taps.emplace_back(p);
        }
        else if(p->Detector & Detector_t::Any_t::CB_Apparatus) {
            cands_cb.emplace_back(p);
            t.b_CBSumVetoE += p->VetoEnergy;
        }
    }

    t.b_nTAPS = unsigned(cands_taps.size());

    if(t.b_nTAPS == 0)
        return;

    steps->Fill("nTAPS>0",1.0);

    t.b_nCB = unsigned(cands_cb.size());

    t.b_CBAvgTime = triggersimu.GetRefTiming();

    if(!isfinite(t.b_CBAvgTime))
        return;
    steps->Fill("CBAvgTime ok",1.0);



    for(const TTaggerHit& taggerhit : event.Reconstructed().TaggerHits) {
        promptrandom.SetTaggerTime(triggersimu.GetCorrectedTaggerTime(taggerhit));
        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;


        t.b_TaggW = promptrandom.FillWeight();
        t.b_TaggE = taggerhit.PhotonEnergy;
        t.b_TaggT = taggerhit.Time;
        t.b_TaggCh = taggerhit.Channel;


        t.b_FitChi2 = 20.0;
        bool kinFit_ok = false;

        for(const auto& it_proton : cands.get_iter()) {

            const TParticlePtr proton = make_shared<TParticle>(ParticleTypeDatabase::Proton, it_proton);


            auto makePhotons = [] (decltype (it_proton)& proton, const TCandidateList& cands) {
                TParticleList l;
                for(const auto& it_photon : cands.get_iter()) {

                    if(it_photon == proton)
                        continue;
                    l.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::Photon, it_photon));
                }

                return l;
            };

            const TParticleList photons = makePhotons(it_proton, cands);

            if(photons.size()==0 || !multiplicities.Contains(unsigned(photons.size())))
                continue;
            steps->Fill("Multiplicity ok",1.0);

            if(save_events)
                manager.SaveEvent();

            LorentzVec photon_sum;
            for(const auto& p : photons) {
                photon_sum += *p;
            }

            // simple missing mass cut
            const auto beam_target = taggerhit.GetPhotonBeam() + LorentzVec({0, 0, 0}, ParticleTypeDatabase::Proton.Mass());
            const auto missing = beam_target - photon_sum;

            if(!interval<double>(750.0,1200.0).Contains(missing.M()))
                continue;

            // proton coplanarity
            const auto copl = std_ext::radian_to_degree(vec2::Phi_mpi_pi(proton->Phi() - photon_sum.Phi() - M_PI ));

            // find the taggerhit with the best E-p conservation Chi2
            utils::KinFitter& fitter = *fitters.at(photons.size()-1);

            // do kinfit
            auto fit_result = fitter.DoFit(taggerhit.PhotonEnergy, proton, photons);


            if(fit_result.Status == APLCON::Result_Status_t::Success
               && fit_result.ChiSquare < t.b_FitChi2) {
                kinFit_ok = true;


                t.b_FitChi2 = fit_result.ChiSquare;
                t.b_NFitIterations = unsigned(fit_result.NIterations);

                t.b_FittedTaggE = fitter.GetFittedBeamE();

                auto fitted_photons = fitter.GetFittedPhotons();
                auto fitted_proton = fitter.GetFittedProton();

                t.b_FittedProton = *fitted_proton;

                t.b_PhotonSum = photon_sum;
                t.b_Missing = missing;

                t.photon_time().clear();
                t.photon_tof().clear();
                t.photon_beta().clear();
                t.photon_PSA().clear();
                t.FittedPhotons().clear();

                t.b_FittedPhotonSum().SetPxPyPzE(0,0,0,0);
                for(const auto& p : fitted_photons) {

                    t.b_FittedPhotonSum() += *p;

                    t.FittedPhotons().push_back(*p);

                    const auto& cluster = p->Candidate->FindCaloCluster();

                    t.photon_PSA().push_back(TVector2(cluster->ShortEnergy, cluster->Energy));
                    t.photon_time().push_back(cluster->Time);

                    if(p->Candidate->Detector & Detector_t::Type_t::TAPS) {

                        const auto dt = taps_detector->GetTimeOfFlight(cluster->Time,
                                                                       cluster->CentralElement,
                                                                       triggersimu.GetRefTiming());
                        t.photon_tof().push_back(dt);

                        const auto beta = taps_detector->GetBeta(*(p->Candidate), triggersimu.GetRefTiming());
                        t.photon_beta().push_back(beta);

                    }
                    else {
                        t.photon_beta().push_back(std_ext::NaN);
                        t.photon_tof().push_back(std_ext::NaN);
                    }

                } // photon loop

                t.b_Proton() = *proton;

                t.b_Proton_vetoE = proton->Candidate->VetoEnergy;

                if(proton->Candidate->Detector & Detector_t::Type_t::TAPS) {

                    const auto taps_cluster = proton->Candidate->FindCaloCluster();
                    const auto dt = taps_detector->GetTimeOfFlight(taps_cluster->Time,
                                                                   taps_cluster->CentralElement,
                                                                   triggersimu.GetRefTiming());
                    const auto beta = taps_detector->GetBeta(*(proton->Candidate), triggersimu.GetRefTiming());


                    t.b_ProtonBeta = beta;
                    t.b_ProtonToF  = dt;
                    t.ProtonPSA() = TVector2(taps_cluster->ShortEnergy, taps_cluster->Energy);

                }

                t.b_ProtonCopl = copl;
                t.b_FittedProtonCopl = std_ext::radian_to_degree(vec2::Phi_mpi_pi(t.b_FittedProton().Phi() - t.b_FittedPhotonSum().Phi() - M_PI ));
            }

            if(kinFit_ok)
                tree->Fill();

        } //proton loop

    } //tagger hit

}

void JustParticles::ShowResult()
{
    tree->Draw("FittedPhotonSum.M()","FitChi2<300 && nTAPS+nCB==3");
}

AUTO_REGISTER_PHYSICS(JustParticles)
