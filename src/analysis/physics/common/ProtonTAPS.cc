#include "ProtonTAPS.h"
#include "base/std_ext/math.h"

#include "TTree.h"

using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace std;


ProtonTAPS::ProtonTAPS(const string& name, OptionsPtr opts):
    Physics(name, opts)
{

    tree = HistFac.makeTTree("tree");


    tree->Branch("nCB",    &b_nCB);
    tree->Branch("nTAPS",  &b_nTAPS);
    tree->Branch("CBAvgTime", &b_CBAvgTime);
    tree->Branch("CBAvgVetoE", &b_CBAvgVetoE);

    tree->Branch("Proton", &b_Proton);

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

void ProtonTAPS::ProcessEvent(const TEvent& event, manager_t&)
{

    TCandidateList cands_taps;
    TCandidateList cands_cb;

    b_CBAvgVetoE = 0;
    for(const auto& p : event.Reconstructed->Candidates) {
        if(p->Detector & Detector_t::Any_t::TAPS_Apparatus) {
            cands_taps.emplace_back(p);
        }
        else if(p->Detector & Detector_t::Any_t::CB_Apparatus) {
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

    // find the proton in TAPS
    b_Proton.Time = numeric_limits<double>::quiet_NaN();
    for(const TCandidatePtr& cand_proton : cands_taps) {
        if(!isfinite(b_Proton.Time) || b_Proton.Time < cand_proton->Time) {
            b_Proton = *cand_proton;
        }
    }

    tree->Fill();
}

void ProtonTAPS::ShowResult()
{
    tree->Draw("Proton.CaloEnergy:(Proton.Time-CBAvgTime)","nCB>2 && CBAvgVetoE < 0.25","colz");
}

AUTO_REGISTER_PHYSICS(ProtonTAPS)
