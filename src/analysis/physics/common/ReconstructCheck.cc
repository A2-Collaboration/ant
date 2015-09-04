#include "physics/common/ReconstructCheck.h"
#include "plot/root_draw.h"
#include "data/Particle.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include <cmath>
#include <iostream>

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace ant::analysis::data;

ReconstructCheck::ReconstructCheck(PhysOptPtr opts):
    Physics("ReconstructCheck",opts),
    nPerEvent(HistFac,"")
{
    const BinSettings e(max(1000.0, atof(GetOption("Emax").c_str())));
    EnergyRec_cb = HistFac.makeTH2D("Energry Reconstruction CB","E_{true} [MeV]","E_{rec} [MeV]", e, e, "Energy_rec_cb");
    EnergyRec_taps = HistFac.makeTH2D("Energry Reconstruction TAPS","E_{true} [MeV]","E_{rec} [MeV]", e, e, "Energy_rec_taps");
}

void ReconstructCheck::ProcessEvent(const Event &event)
{
    if(event.MCTrue().Particles().GetAll().size() == 1) {
        const auto& p = event.MCTrue().Particles().GetAll().at(0);

        if( event.Reconstructed().Candidates().size() == 1)
        for(auto& cand : event.Reconstructed().Candidates()) {

            if(cand->Detector() & Detector_t::Any_t::CB) {
                EnergyRec_cb->Fill(p->Ek(), cand->ClusterEnergy());
            } else if( cand->Detector() & Detector_t::Any_t::TAPS) {
                EnergyRec_taps->Fill(p->Ek(), cand->ClusterEnergy());
            }

        }


        nPerEvent.Fill(p, event.Reconstructed().Candidates());

        auto entry = nPerEvent_type.find(&(p->Type()));
        candidatesEvent_t* c = nullptr;
        if(entry == nPerEvent_type.end()) {
            auto pos = nPerEvent_type.insert(make_pair( &(p->Type()), candidatesEvent_t(HistFac,p->Type().Name())));
            c = &(pos.first->second);
        } else {
            c = &(entry->second);
        }

        c->Fill(p,event.Reconstructed().Candidates());

    }
}


void ReconstructCheck::Finish()
{

}

canvas& draw(canvas& c, const std::list<TH1*>& l) {
    for(auto& h : l) {
        c << h;
    }
    return c;
}

void ReconstructCheck::ShowResult()
{
    canvas("ReconstructCheck")
            << drawoption("colz")
            << EnergyRec_cb << EnergyRec_taps << endc;

    canvas candidates("Candidates");

    candidates << drawoption("colz");

    draw(candidates , nPerEvent.Hists());

    if(nPerEvent_type.size() > 1) {
        for(auto& e : nPerEvent_type) {
            draw(candidates , e.second.Hists());
        }
    }

    candidates << endc;
}




ReconstructCheck::candidatesEvent_t::candidatesEvent_t(SmartHistFactory& f, const string& prefix)
{
    nPerEvent     = f.makeTH1D(prefix+" Candidates/Event", "Candidates/Event","",BinSettings(15),prefix+"candEvent");
    nPerEventPerE = f.makeTH2D(prefix+" Candidates/Event/Energy","MC True Energy [MeV]","Candidates/Event",BinSettings(1000),BinSettings(15),prefix+"candEventEnergy");
    splitPerEvent = f.makeTH1D(prefix+" Split Clusters/Event", "# split clusters/Event","",BinSettings(30),prefix+"splits");
    splitPos      = f.makeTH2D(prefix+" Pos of split clusters","cos(#theta)","#phi",BinSettings(360,-1,1),BinSettings(360,-180,180),prefix+"splitpos");
}

void ReconstructCheck::candidatesEvent_t::Fill(const ParticlePtr& mctrue, const CandidateList& cand)
{
    nPerEvent->Fill(cand.size());
    nPerEventPerE->Fill(mctrue->Ek(), cand.size());

    unsigned nsplit(0);

    for(const CandidatePtr& c : cand) {
        for(const Cluster& cl : c->Clusters) {
            if(cl.flags.isChecked(Cluster::Flag::Split)) {
                ++nsplit;
                splitPos->Fill(cos(mctrue->Theta()), mctrue->Phi()*TMath::RadToDeg());
            }
        }
    }
    splitPerEvent->Fill(nsplit);
}

std::list<TH1*> ReconstructCheck::candidatesEvent_t::Hists()
{
    return {nPerEvent, nPerEventPerE, splitPerEvent, splitPos};
}

AUTO_REGISTER_PHYSICS(ReconstructCheck, "ReconstructCheck")
