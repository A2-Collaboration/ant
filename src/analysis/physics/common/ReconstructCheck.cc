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

    canvas splitanglesc("Split Angles");
    hstack s("splitanglesstack");
    for(auto h : nPerEvent.mult2_split_angles) {
        s << (TH1D*)h;
    }

    splitanglesc << drawoption("nostack") << s << endc;

    canvas c_pid("PID");
    c_pid << nPerEvent.nCharged << drawoption("colz") << nPerEvent.cluserSize_CB << nPerEvent.cluserSize_TAPS
          << nPerEvent.dEE_CB << nPerEvent.dEE_CB_true
          << nPerEvent.dEE_TAPS << nPerEvent.dEE_TAPS_true
          << nPerEvent.posCharged
          << endc;
}

void LabelBins(TAxis* x) {

    for(Int_t i=1; i<=x->GetNbins(); ++i) {
        const auto a = x->GetBinLowEdge(i);
        x->SetBinLabel(i,to_string(int(a)).c_str());
    }
}


ReconstructCheck::candidatesEvent_t::candidatesEvent_t(SmartHistFactory& f, const string& prefix)
{
    const BinSettings energy(1000);
    const BinSettings vetoEnergy(100,0,25);
    const BinSettings clusersize(19,1,20);
    const BinSettings costheta(360,-1,1);
    const BinSettings phi(360,-180,180);

    nPerEvent     = f.makeTH1D(prefix+" Candidates/Event", "Candidates/Event","",BinSettings(10),prefix+"candEvent");
    LabelBins(nPerEvent->GetXaxis());
    nPerEvent->SetFillColor(kGray);

    nPerEventPerE = f.makeTH2D(prefix+" Candidates/Event/Energy","MC True Energy [MeV]","Candidates/Event",BinSettings(1000),BinSettings(10),prefix+"candEventEnergy");
    LabelBins(nPerEventPerE->GetYaxis());

    splitPerEvent = f.makeTH1D(prefix+" Split-flagged Clusters/Event", "# split clusters/Event","",BinSettings(20),prefix+"splits");
    LabelBins(splitPerEvent->GetXaxis());
    splitPerEvent->SetFillColor(kGray);

    splitPos      = f.makeTH2D(prefix+" Pos of split-flagged clusters","cos(#theta)","#phi",costheta,phi,prefix+"splitpos");
    splitPos->SetStats(false);

    multiplicity_map      = f.makeTH3D(prefix+" Multitplicity Map","cos(#theta_{True})","#phi_{True}","Mult",costheta,phi,BinSettings(10),prefix+"mults");
    multiplicity_map->SetStats(false);

    mult2_split_angles.resize(3);
    for(size_t i=0;i<mult2_split_angles.size();++i) {
        mult2_split_angles[i] = f.makeTH1D(prefix+"Mult==2 cluster angle "+to_string(i),"#alpha [#circ]","",BinSettings(180,0,90),prefix+"mult2_"+to_string(i));
    }

    cluserSize_TAPS = f.makeTH2D(prefix+" Cluster Size TAPS","E_{True} [MeV]","Elements",energy, clusersize,prefix+"_clustersize_TAPS");
    LabelBins(cluserSize_TAPS->GetYaxis());

    cluserSize_CB   = f.makeTH2D(prefix+" Cluster Size CB",  "E_{True} [MeV]","Elements",energy, clusersize,prefix+"_clustersize_CB");
    LabelBins(cluserSize_CB->GetYaxis());


    dEE_TAPS = f.makeTH2D(prefix+" dEE TAPS", "E [MeV]","VetoEnergy [MeV]",energy, vetoEnergy,prefix+"_dEE_TAPS");
    dEE_TAPS_true = f.makeTH2D(prefix+" dEE TAPS (true E)", "E_{True} [MeV]","VetoEnergy [MeV]",energy, vetoEnergy,prefix+"_dEE_true_TAPS");

    dEE_CB = f.makeTH2D(prefix+" dEE CB", "E [MeV]","VetoEnergy [MeV]",energy, vetoEnergy,prefix+"_dEE_CB");
    dEE_CB_true = f.makeTH2D(prefix+" dEE CB (True E)", "E_{True} [MeV]","VetoEnergy [MeV]",energy, vetoEnergy,prefix+"_dEE_true_CB");

    nCharged        = f.makeTH1D(prefix+" N Charged (VetoEnergy > 0)", "# charged candidates", "", BinSettings(10),prefix+"_ncharged");
    LabelBins(nCharged->GetXaxis());
    nCharged->SetFillColor(kGray);

    posCharged = f.makeTH2D(prefix+"Position of charged Candidates","cos(#theta_{True})","#phi [#circ]",costheta,phi,prefix+"_posCharged");

}

double angle(const data::Candidate& c1, const data::Candidate& c2) {
    TVector3 v1(1,0,0);
    v1.SetTheta(c1.Theta());
    v1.SetPhi(c1.Phi());

    TVector3 v2(1,0,0);
    v2.SetTheta(c2.Theta());
    v2.SetPhi(c2.Phi());

    return v1.Angle(v2);
}

void ReconstructCheck::candidatesEvent_t::Fill(const ParticlePtr& mctrue, const CandidateList& cand)
{
    const auto mc_phi = mctrue->Phi()*TMath::RadToDeg();
    const auto mc_cos_theta = cos(mctrue->Theta());
    const auto mc_energy = mctrue->Ek();

    nPerEvent->Fill(cand.size());
    nPerEventPerE->Fill(mctrue->Ek(), cand.size());

    unsigned nsplit(0);
    unsigned ncharged(0);

    for(const CandidatePtr& c : cand) {
        if(c->VetoEnergy() > 0.0) {
            ncharged++;
            posCharged->Fill(mc_cos_theta,mc_phi);
        }

        if(c->Detector() & Detector_t::Any_t::CB) {
            cluserSize_CB->Fill(mc_energy, c->ClusterSize());
            dEE_CB->Fill(c->ClusterEnergy(), c->VetoEnergy());
            dEE_CB_true->Fill(mc_energy, c->VetoEnergy());
        } else if(c->Detector() & Detector_t::Any_t::TAPS) {
            cluserSize_TAPS->Fill(mc_energy, c->ClusterSize());
            dEE_TAPS->Fill(c->ClusterEnergy(), c->VetoEnergy());
            dEE_TAPS_true->Fill(mc_energy, c->VetoEnergy());
        }

        for(const Cluster& cl : c->Clusters) {
            if(cl.flags.isChecked(Cluster::Flag::Split)) {
                ++nsplit;
                splitPos->Fill(mc_cos_theta, mc_phi);
            }
        }
    }
    nCharged->Fill(ncharged);
    splitPerEvent->Fill(nsplit);
    multiplicity_map->Fill(mc_cos_theta, mc_phi, cand.size());

    if(cand.size() == 2) {
        if(nsplit<3) {
            mult2_split_angles[nsplit]->Fill(angle(*cand.at(0), *cand.at(1))*TMath::RadToDeg());
        }
    }

}

std::list<TH1*> ReconstructCheck::candidatesEvent_t::Hists()
{
    return {nPerEvent, nPerEventPerE, splitPerEvent, splitPos, multiplicity_map};
}

AUTO_REGISTER_PHYSICS(ReconstructCheck, "ReconstructCheck")
