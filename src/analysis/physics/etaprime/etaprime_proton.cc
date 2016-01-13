#include "etaprime_proton.h"
#include "base/std_ext/math.h"

#include "TTree.h"

using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::data;
using namespace ant::analysis::physics;
using namespace std;


EtapProton::EtapProton(const string& name, PhysOptPtr opts):
    Physics(name, opts)
{

    tree = HistFac.makeTTree("tree");


    tree->Branch("nCB",    &b_nCB);
    tree->Branch("nTAPS",  &b_nTAPS);
    tree->Branch("CBAvgTime", &b_CBAvgTime);
    tree->Branch("CBAvgVetoE", &b_CBAvgVetoE);

    tree->Branch("Proton", &b_Proton);

    tree->Branch("PhotonSum", &b_PhotonSum);
    tree->Branch("ProtonCopl", &b_ProtonCopl);

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

    data::CandidateList cands_taps;
    data::CandidateList cands_cb;

    b_CBAvgVetoE = 0;
    for(const auto& p : event.Reconstructed.Candidates) {
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
    b_CBAvgVetoE /= b_nCB;

    // find the proton candidate in TAPS, ie. the slowest cluster in TAPS
    /// \todo this could be improved with some proper beta calculation,
    /// but that needs precise TAPS timing offsets for each element:
    ///
    /// beta = s / (s + c*dt*cos(theta))
    ///
    /// dt = 0 if particle is photon (that's crucial!)
    b_Proton.Time = numeric_limits<double>::quiet_NaN();
    data::CandidatePtr proton;
    for(const CandidatePtr& cand_taps : cands_taps) {
        if(!isfinite(b_Proton.Time) || b_Proton.Time < cand_taps->Time) {
            b_Proton = *cand_taps;
            proton = cand_taps;
        }
    }

    // create "photons" from all other clusters
    vector<Particle> photons;
    for(const auto& cand_cb : cands_cb)
        photons.emplace_back(ParticleTypeDatabase::Photon, cand_cb);
    for(const auto& cand_taps : cands_taps) {
        if(cand_taps != proton)
            photons.emplace_back(ParticleTypeDatabase::Photon, cand_taps);
    }

    b_PhotonSum.SetPxPyPzE(0,0,0,0);
    for(const auto& p : photons) {
       b_PhotonSum += p;
    }

    // proton coplanarity
    b_ProtonCopl = std_ext::radian_to_degree(TVector2::Phi_mpi_pi(proton->Phi - b_PhotonSum.Phi() - M_PI ));

    tree->Fill();
}

void EtapProton::ShowResult()
{
}

AUTO_REGISTER_PHYSICS(EtapProton)
