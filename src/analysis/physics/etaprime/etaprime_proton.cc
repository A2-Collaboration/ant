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
    promptrandom.AddPromptRange({-2.5,1.5}); // slight offset due to CBAvgTime reference
    promptrandom.AddRandomRange({-50,-10});  // just ensure to be way off prompt peak
    promptrandom.AddRandomRange({  10,50});

    tree = HistFac.makeTTree("tree");


    tree->Branch("nCB",    &b_nCB);
    tree->Branch("nTAPS",  &b_nTAPS);
    tree->Branch("CBAvgTime", &b_CBAvgTime);
    tree->Branch("CBAvgVetoE", &b_CBAvgVetoE);

    tree->Branch("Proton", &b_Proton);
    tree->Branch("Photons", &b_Photons);
    tree->Branch("PhotonSum", &b_PhotonSum);
    tree->Branch("ProtonCopl", &b_ProtonCopl);

    tree->Branch("BestChi2", &b_BestChi2);
    tree->Branch("TaggW",  &b_TaggW);
    tree->Branch("TaggE",  &b_TaggE);
    tree->Branch("TaggT",  &b_TaggT);
    tree->Branch("TaggCh", &b_TaggCh);

    tree->Branch("FittedProton", &b_FittedProton);
    tree->Branch("FittedPhotons", &b_FittedPhotons);
    tree->Branch("FittedPhotonSum", &b_FittedPhotonSum);
    tree->Branch("FittedProtonCopl", &b_FittedProtonCopl);

    // prepare fitters for all multiplicities
    const auto setup = ant::ExpConfig::Setup::GetLastFound();
    for(unsigned mult=1;mult<=11;mult++) {
        auto fitter = std_ext::make_unique<utils::KinFitter>(
                          std_ext::formatter() << GetName() << mult,
                          mult
                          );
        fitter->LoadSigmaData(setup->GetPhysicsFilesDirectory()+"/FitterSigmas.root");
        fitters.emplace_back(move(fitter));
    }
}

template <typename T>
double TimeAverage(const T& cands) {
    double time   = 0.0;
    double energy = 0.0;
    for(const auto& c : cands) {
        time += c->Time * c->CaloEnergy;
        energy += c->CaloEnergy;
    }
    return time / energy;
}

void EtapProton::ProcessEvent(const data::Event& event)
{

    const auto& cands = event.Reconstructed.Candidates;
    data::CandidateList cands_taps;
    data::CandidateList cands_cb;

    b_CBAvgVetoE = 0;
    for(const auto& p : cands) {
        if(p->GetDetector() & Detector_t::Any_t::TAPS_Apparatus) {
            cands_taps.emplace_back(p);
        }
        else if(p->GetDetector() & Detector_t::Any_t::CB_Apparatus) {
            cands_cb.emplace_back(p);
            b_CBAvgVetoE += p->VetoEnergy;
        }
    }
    b_nTAPS = cands_taps.size();
    if(b_nTAPS == 0)
        return;

    b_nCB = cands_cb.size();
    b_CBAvgTime = TimeAverage(cands_cb);
    if(!isfinite(b_CBAvgTime))
        return;
    b_CBAvgVetoE /= b_nCB;

    // find the proton candidate in TAPS, ie. the slowest cluster in TAPS
    /// \todo this could be improved with some proper beta calculation,
    /// but that needs precise TAPS timing offsets for each element:
    ///
    /// beta = s / (s + c*dt*cos(theta))
    ///
    /// dt = 0 if particle is photon (that's crucial!)
    b_Proton.Time = numeric_limits<double>::quiet_NaN();
    data::ParticlePtr proton;
    for(const CandidatePtr& cand_taps : cands_taps) {
        if(!isfinite(b_Proton.Time) || b_Proton.Time < cand_taps->Time) {
            b_Proton = *cand_taps;
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

    if(photons.size()==0)
        return;

    b_PhotonSum.SetPxPyPzE(0,0,0,0);
    b_Photons.resize(0);
    for(const auto& p : photons) {
       b_PhotonSum += *p;
       b_Photons.emplace_back(*p->Candidate);
    }

    // proton coplanarity
    b_ProtonCopl = std_ext::radian_to_degree(TVector2::Phi_mpi_pi(proton->Phi() - b_PhotonSum.Phi() - M_PI ));

    // find the taggerhit with the best E-p conservation Chi2
    if(photons.size()>fitters.size())
        return;
    auto& fitter = *fitters[photons.size()-1];
    b_BestChi2 = std::numeric_limits<double>::quiet_NaN();
    for(const data::TaggerHit& taggerhit : event.Reconstructed.TaggerHits) {
        promptrandom.SetTaggerHit(taggerhit.Time - b_CBAvgTime);
        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;

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

        // only update stuff if better chi2 found
        if(!isfinite(b_BestChi2) || b_BestChi2 > fit_result.ChiSquare) {
            b_BestChi2 = fit_result.ChiSquare;

            b_TaggW = promptrandom.FillWeight();
            b_TaggE = taggerhit.PhotonEnergy;
            b_TaggT = taggerhit.Time;
            b_TaggCh = taggerhit.Channel;

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

    // ignore events without successful fit
    if(isfinite(b_BestChi2))
        tree->Fill();
}

void EtapProton::ShowResult()
{
    tree->Draw("FittedPhotonSum.M()","BestChi2<300 && CBAvgVetoE<0.25 && @Photons.size()==2");
}

AUTO_REGISTER_PHYSICS(EtapProton)
