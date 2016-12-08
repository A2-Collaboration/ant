#include "JustParticles.h"
#include "base/std_ext/math.h"
#include "expconfig/ExpConfig.h"

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

    promptrandom.AddPromptRange({-2.5,1.5}); // slight offset due to CBAvgTime reference
    promptrandom.AddRandomRange({-50,-10});  // just ensure to be way off prompt peak
    promptrandom.AddRandomRange({  10,50});

    tree = HistFac.makeTTree("tree");

    // prepare fitters for all multiplicities
    fitters.resize(enclosing.Stop());
    const auto setup = ant::ExpConfig::Setup::GetLastFound();
    if(!setup)
        throw runtime_error("EtapProton needs a setup");

    auto uncertainty = make_shared<utils::UncertaintyModels::FitterSergey>();

    for(unsigned mult=enclosing.Start();mult<=enclosing.Stop();mult++) {
        if(!multiplicities.Contains(mult))
            continue;
        auto fitter = std_ext::make_unique<utils::KinFitter>(
                          std_ext::formatter() << "FitMult" << mult,
                          mult,
                          uncertainty
                          );

        fitters[mult-1] = move(fitter);
    }

    // needed for calculating ToF
    taps_detector = ExpConfig::Setup::GetDetector<expconfig::detector::TAPS>();
    if(!taps_detector)
        throw runtime_error("EtapProton needs TAPS detector in setup");
}

void JustParticles::ProcessEvent(const TEvent& event, manager_t& manager)
{
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
    t.b_CBAvgTime = event.Reconstructed().Trigger.CBTiming;
    if(!isfinite(t.b_CBAvgTime))
        return;
    steps->Fill("CBAvgTime ok",1.0);


    // find the proton candidate in TAPS, ie. the lowest beta=v/c in TAPS
    t.b_ProtonBeta = numeric_limits<double>::quiet_NaN();
    TParticlePtr proton;
    for(const TCandidatePtr& cand_taps : cands_taps) {
        // calculate the beta = v/c of the particle from time of flight
        // note that the time of flight is only correct if the correct reference time
        // is used...
        const auto& trigger_reftime = event.Reconstructed().Trigger.CBTiming;
        const auto taps_cluster = cand_taps->FindCaloCluster();
        const auto dt = taps_detector->GetTimeOfFlight(taps_cluster->Time,
                                                       taps_cluster->CentralElement,
                                                       trigger_reftime);
        const auto beta = taps_detector->GetBeta(*cand_taps, trigger_reftime);

        if(!isfinite(t.b_ProtonBeta) || t.b_ProtonBeta > beta) {
            t.b_ProtonBeta = beta;
            t.b_ProtonToF  = dt;
            t.b_ProtonPSA_R     = taps_cluster->GetPSARadius();
            t.b_ProtonPSA_Angle = taps_cluster->GetPSAAngle();
            proton = make_shared<TParticle>(ParticleTypeDatabase::Proton, cand_taps);
            t.b_Proton = *proton;
            t.b_Proton_vetoE = cand_taps->VetoEnergy;
        }

    }

    // create "photons" from all other clusters
    TParticleList photons;
    for(const auto& cand_cb : cands_cb)
        photons.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::Photon, cand_cb));
    for(const auto& cand_taps : cands_taps) {
        if(cand_taps != proton->Candidate)
            photons.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::Photon, cand_taps));
    }

    assert(photons.size() == cands.size()-1);

    if(photons.size()==0 || !multiplicities.Contains(unsigned(photons.size())))
        return;
    steps->Fill("Multiplicity ok",1.0);

    if(save_events)
        manager.SaveEvent();

    t.b_PhotonSum().SetPxPyPzE(0,0,0,0);
    for(const auto& p : photons) {
       t.b_PhotonSum() += *p;
    }

    // proton coplanarity
    t.b_ProtonCopl = std_ext::radian_to_degree(vec2::Phi_mpi_pi(proton->Phi() - t.b_PhotonSum().Phi() - M_PI ));

    // find the taggerhit with the best E-p conservation Chi2
    utils::KinFitter& fitter = *fitters.at(photons.size()-1);

    bool kinFit_ok = false;
    for(const TTaggerHit& taggerhit : event.Reconstructed().TaggerHits) {
        promptrandom.SetTaggerHit(taggerhit.Time - t.b_CBAvgTime);
        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        // simple missing mass cut
        const auto beam_target = taggerhit.GetPhotonBeam() + LorentzVec({0, 0, 0}, ParticleTypeDatabase::Proton.Mass());
        t.b_Missing = beam_target - t.b_PhotonSum();

        t.b_TaggW = promptrandom.FillWeight();
        t.b_TaggE = taggerhit.PhotonEnergy;
        t.b_TaggT = taggerhit.Time;
        t.b_TaggCh = taggerhit.Channel;

        // do kinfit
        fitter.SetEgammaBeam(taggerhit.PhotonEnergy);
        fitter.SetProton(proton);
        fitter.SetPhotons(photons);
        auto fit_result = fitter.DoFit();


        t.b_FitStatus = static_cast<unsigned>(fit_result.Status);
        t.b_FitChi2 = numeric_limits<double>::quiet_NaN();
        t.b_NFitIterations = 0;
        t.b_FittedTaggE = numeric_limits<double>::quiet_NaN();
        t.b_FittedPhotonSum().SetPxPyPzE(0,0,0,0);
        t.b_FittedProton().SetPxPyPzE(0,0,0,0);
        t.b_FittedProtonCopl = numeric_limits<double>::quiet_NaN();

        if(fit_result.Status == APLCON::Result_Status_t::Success) {
            kinFit_ok = true;

            t.b_FitChi2 = fit_result.ChiSquare;
            t.b_NFitIterations = unsigned(fit_result.NIterations);

            t.b_FittedTaggE = fitter.GetFittedBeamE();

            auto fitted_photons = fitter.GetFittedPhotons();
            auto fitted_proton = fitter.GetFittedProton();

            t.b_FittedProton = *fitted_proton;

            for(const auto& p : fitted_photons) {
                t.b_FittedPhotonSum() += *p;
            }

            t.b_FittedProtonCopl = std_ext::radian_to_degree(vec2::Phi_mpi_pi(t.b_FittedProton().Phi() - t.b_FittedPhotonSum().Phi() - M_PI ));
        }

        tree->Fill();
    }

    if(kinFit_ok)
        steps->Fill("KinFit ok",1.0);
}

void JustParticles::ShowResult()
{
    tree->Draw("FittedPhotonSum.M()","FitChi2<300 && nTAPS+nCB==3");
}

AUTO_REGISTER_PHYSICS(JustParticles)
