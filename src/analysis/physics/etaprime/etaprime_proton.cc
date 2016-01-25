#include "etaprime_proton.h"
#include "base/std_ext/math.h"
#include "expconfig/ExpConfig.h"

#include "TTree.h"
#include "TLorentzVector.h"

#include <cassert>

using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::data;
using namespace ant::analysis::physics;
using namespace std;


EtapProton::EtapProton(const string& name, PhysOptPtr opts):
    Physics(name, opts)
{
    multiplicities = opts->Get<decltype(multiplicities)>("PhotonMulti",{{2,2},{4,4},{6,6}});
    auto enclosing = multiplicities.EnclosingInterval();
    if(!enclosing.IsSane() || enclosing.Start()<1)
        throw runtime_error("Given photon multiplicities not sane");


    steps = HistFac.makeTH1D("Steps","","",BinSettings(10),"steps");

    promptrandom.AddPromptRange({-2.5,1.5}); // slight offset due to CBAvgTime reference
    promptrandom.AddRandomRange({-50,-10});  // just ensure to be way off prompt peak
    promptrandom.AddRandomRange({  10,50});

    tree = HistFac.makeTTree("tree");

    tree->Branch("nCB",        &b_nCB);
    tree->Branch("nTAPS",      &b_nTAPS);
    tree->Branch("CBAvgTime",  &b_CBAvgTime);
    tree->Branch("CBSumVetoE", &b_CBSumVetoE);

    tree->Branch("Proton",     &b_Proton);
    tree->Branch("Photons",    &b_Photons);
    tree->Branch("PhotonSum",  &b_PhotonSum);
    tree->Branch("ProtonCopl", &b_ProtonCopl);
    tree->Branch("ProtonBeta", &b_ProtonBeta);
    tree->Branch("ProtonToF",  &b_ProtonToF);


    tree->Branch("BestChi2",  &b_BestChi2);
    tree->Branch("NGoodFits", &b_NGoodFits); // number of good fits
    tree->Branch("NFitIterations", &b_NFitIterations); // number of good fits

    tree->Branch("TaggW",  &b_TaggW); // prompt/random weight
    tree->Branch("TaggE",  &b_TaggE);
    tree->Branch("TaggT",  &b_TaggT);
    tree->Branch("TaggCh", &b_TaggCh);
    tree->Branch("TaggN",  &b_TaggN); // number of hits


    tree->Branch("FittedTaggE",      &b_FittedTaggE);
    tree->Branch("FittedProton",     &b_FittedProton);
    tree->Branch("FittedPhotons",    &b_FittedPhotons);
    tree->Branch("FittedPhotonSum",  &b_FittedPhotonSum);
    tree->Branch("FittedProtonCopl", &b_FittedProtonCopl);

    // prepare fitters for all multiplicities
    fitters.resize(enclosing.Stop());
    const auto setup = ant::ExpConfig::Setup::GetLastFound();
    for(unsigned mult=enclosing.Start();mult<=enclosing.Stop();mult++) {
        if(!multiplicities.Contains(mult))
            continue;
        auto fitter = std_ext::make_unique<utils::KinFitter>(
                          std_ext::formatter() << "FitMult" << mult,
                          mult
                          );
        fitter->LoadSigmaData(setup->GetPhysicsFilesDirectory()+"/FitterSigmas.root");
        fitter->SetupBranches(tree);
        fitters[mult-1] = move(fitter);
    }

    // needed for calculating ToF
    taps_detector = ExpConfig::Setup::GetDetector<expconfig::detector::TAPS>();
}

void EtapProton::ProcessEvent(const data::Event& event)
{
    steps->Fill("Seen",1.0);

    const auto& cands = event.Reconstructed.Candidates;
    data::CandidateList cands_taps;
    data::CandidateList cands_cb;

    b_CBSumVetoE = 0;
    for(const auto& p : cands) {
        if(p->GetDetector() & Detector_t::Any_t::TAPS_Apparatus) {
            cands_taps.emplace_back(p);
        }
        else if(p->GetDetector() & Detector_t::Any_t::CB_Apparatus) {
            cands_cb.emplace_back(p);
            b_CBSumVetoE += p->VetoEnergy;
        }
    }
    b_nTAPS = cands_taps.size();
    if(b_nTAPS == 0)
        return;
    steps->Fill("nTAPS>0",1.0);

    b_nCB = cands_cb.size();
    b_CBAvgTime = event.Reconstructed.Trigger.CBTiming;
    if(!isfinite(b_CBAvgTime))
        return;
    steps->Fill("CBAvgTime ok",1.0);


    // find the proton candidate in TAPS, ie. the lowest beta=v/c in TAPS
    b_ProtonBeta = numeric_limits<double>::quiet_NaN();
    data::ParticlePtr proton;
    for(const CandidatePtr& cand_taps : cands_taps) {
        // calculate the beta = v/c of the particle from time of flight
        // note that the time of flight is only correct if the correct reference time
        // is used...
        auto taps_cluster = cand_taps->FindCaloCluster();
        const double dt = taps_detector->GetTimeOfFlight(taps_cluster->Time, taps_cluster->CentralElement,
                                                          event.Reconstructed.Trigger.CBTiming);
        const double s = taps_detector->GetZPosition();
        constexpr double c = 30; // velocity of light in cm/ns

        const double beta = s / (s + c * dt * cos(cand_taps->Theta));

        if(!isfinite(b_ProtonBeta) || b_ProtonBeta > beta) {
            b_Proton = *cand_taps;
            b_ProtonBeta = beta;
            b_ProtonToF = dt;
            proton = make_shared<Particle>(ParticleTypeDatabase::Proton, cand_taps);
        }
    }

    // create "photons" from all other clusters
    data::ParticleList photons;
    for(const auto& cand_cb : cands_cb)
        photons.emplace_back(make_shared<Particle>(ParticleTypeDatabase::Photon, cand_cb));
    for(const auto& cand_taps : cands_taps) {
        if(cand_taps != proton->Candidate)
            photons.emplace_back(make_shared<Particle>(ParticleTypeDatabase::Photon, cand_taps));
    }

    assert(photons.size() == cands.size()-1);

    if(photons.size()==0 || !multiplicities.Contains(photons.size()))
        return;
    steps->Fill("Multiplicity ok",1.0);

    b_PhotonSum.SetPxPyPzE(0,0,0,0);
    b_Photons.resize(0);
    for(const auto& p : photons) {
       b_PhotonSum += *p;
       b_Photons.emplace_back(*p->Candidate);
    }

    // proton coplanarity
    b_ProtonCopl = std_ext::radian_to_degree(TVector2::Phi_mpi_pi(proton->Phi() - b_PhotonSum.Phi() - M_PI ));

    // find the taggerhit with the best E-p conservation Chi2
    utils::KinFitter& fitter = *fitters.at(photons.size()-1);
    b_BestChi2 = std::numeric_limits<double>::quiet_NaN();
    b_TaggW = 0; // ROOT does not work properly if weights are NaN
    b_TaggT = std::numeric_limits<double>::quiet_NaN();
    b_TaggE = std::numeric_limits<double>::quiet_NaN();
    b_FittedTaggE = std::numeric_limits<double>::quiet_NaN();
    b_TaggN = 0;
    b_NGoodFits = 0;
    b_NFitIterations = 0;
    b_FittedProtonCopl = std::numeric_limits<double>::quiet_NaN();
    for(const data::TaggerHit& taggerhit : event.Reconstructed.TaggerHits) {
        promptrandom.SetTaggerHit(taggerhit.Time - b_CBAvgTime);
        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;
        b_TaggN++;

        // simple missing mass cut
        const TLorentzVector beam_target = taggerhit.GetPhotonBeam() + TLorentzVector(0, 0, 0, ParticleTypeDatabase::Proton.Mass());
        const TLorentzVector missing = beam_target - b_PhotonSum;

        // skip tagger hits which don't match found proton
        const double angle_p_calcp = std_ext::radian_to_degree(missing.Angle(proton->Vect()));
        if(angle_p_calcp > 15.0)
            continue;

        // do kinfit
        fitter.SetEgammaBeam(taggerhit.PhotonEnergy);
        fitter.SetProton(proton);
        fitter.SetPhotons(photons);
        auto fit_result = fitter.DoFit();

        if(fit_result.Status != APLCON::Result_Status_t::Success)
            continue;
        b_NGoodFits++;


        // only update stuff if better chi2 found
        if(!isfinite(b_BestChi2) || b_BestChi2 > fit_result.ChiSquare) {
            b_BestChi2 = fit_result.ChiSquare;
            b_NFitIterations = fit_result.NIterations;

            b_TaggW = promptrandom.FillWeight();
            b_TaggE = taggerhit.PhotonEnergy;
            b_TaggT = taggerhit.Time;
            b_TaggCh = taggerhit.Channel;

            b_FittedTaggE = fitter.GetFittedBeamE();

            auto fitted_photons = fitter.GetFittedPhotons();
            auto fitted_proton = fitter.GetFittedProton();

            b_FittedProton = *fitted_proton;

            b_FittedPhotonSum.SetPxPyPzE(0,0,0,0);
            b_FittedPhotons.resize(0);
            for(const auto& p : fitted_photons) {
               b_FittedPhotonSum += *p;
               b_FittedPhotons.emplace_back(*p);
            }

            b_FittedProtonCopl = std_ext::radian_to_degree(TVector2::Phi_mpi_pi(fitted_proton->Phi() - b_FittedPhotonSum.Phi() - M_PI ));
        }

    }

    if(isfinite(b_BestChi2))
        steps->Fill("KinFit ok",1.0);

    tree->Fill();
}

void EtapProton::ShowResult()
{
    tree->Draw("FittedPhotonSum.M()","BestChi2<300 && @Photons.size()==2");
}

AUTO_REGISTER_PHYSICS(EtapProton)
