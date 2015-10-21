#include "physics/common/CandidatesAnalysis.h"
#include "plot/root_draw.h"
#include "data/Particle.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include <cmath>

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace ant::analysis::data;

CandidatesAnalysis::CandidatesAnalysis(const std::string& name, PhysOptPtr opts):
    Physics(name,opts)
{
    nCandidatesEvent = HistFac.makeTH1D("Candidates/Event","Candidates","",BinSettings(15),"nCand");
    CandMultiplicities = HistFac.makeTH1D("Candidates Multi","","",BinSettings(15),"CandMult");

    for(int i=1; i<CandMultiplicities->GetNbinsX(); ++i) {
        CandMultiplicities->GetXaxis()->SetBinLabel(i,string(to_string(i-1)+"+").c_str());
    }

    energy = HistFac.makeTH1D("Energy","E [MeV]","",BinSettings(1000),"energy");
    theta  = HistFac.makeTH1D("Theta","#theta [#circ]","",BinSettings(360,0,180),"theta");
    phi    = HistFac.makeTH1D("Phi","#phi [#circ]","",BinSettings(720,-180,180),"phi");
    ggIM         = HistFac.makeTH1D("2 Neutral Candidates IM","M [MeV]","",BinSettings(1000),"ggIM");
    ttIM         = HistFac.makeTH1D("2 Candidates IM","M [MeV]","",BinSettings(1000),"ttIM");
    cbdEE = HistFac.makeTH2D("CB dE-E","E_{CB} [MeV]","dE_{PID} [MeV]", BinSettings(1000),BinSettings(100,0,30),"cb_dEE");
    cbtof = HistFac.makeTH2D("CB ToF","t_{CB} [ns]","E_{CB} [MeV]", BinSettings(300,-15,15),BinSettings(1000),"cb_tof");
    tapsdEE = HistFac.makeTH2D("TAPS dE-E","E_{TAPS} [MeV]","dE_{TAPSVeto} [MeV]", BinSettings(1000),BinSettings(100,0,30),"taps_dEE");
    tapstof = HistFac.makeTH2D("TAPS ToF","t_{TAPS} [ns]","E_{TAPS} [MeV]", BinSettings(300,-15,15),BinSettings(1000),"taps_tof");
    detectors = HistFac.makeTH1D("Detectors","","", BinSettings(1),"detectors");
}

void CandidatesAnalysis::ProcessEvent(const Event &event)
{
    const auto& candidates = event.Reconstructed().Candidates();

    nCandidatesEvent->Fill(candidates.size());

    for(size_t i=0; i<=candidates.size(); ++i) {
        CandMultiplicities->Fill(i);
    }

    CandidateList::const_iterator i = candidates.begin();

    while(i!=candidates.end()) {
        const CandidatePtr& ci = *i;
        energy->Fill((*i)->ClusterEnergy());
        theta->Fill((*i)->Theta()*TMath::RadToDeg());
        phi->Fill((*i)->Phi()*TMath::RadToDeg());
        detectors->Fill(string(ci->Detector()).c_str(),1);

        if((*i)->ClusterEnergy()>20.0) {

            if(ci->Detector() & Detector_t::Any_t::CB) {
                cbdEE->Fill(ci->ClusterEnergy(),ci->VetoEnergy());
                cbtof->Fill(ci->Time(), ci->ClusterEnergy());
            } else if(ci->Detector() & Detector_t::Any_t::TAPS) {
                tapsdEE->Fill(ci->ClusterEnergy(),ci->VetoEnergy());
                tapstof->Fill(ci->Time(), ci->ClusterEnergy());
            }

            const Particle a(ParticleTypeDatabase::Photon,*i);
            CandidateList::const_iterator j = i;
            ++j;
            while(j!=candidates.end()) {
                if((*i)->ClusterEnergy()>20.0) {
                    const Particle b(ParticleTypeDatabase::Photon,*j);
                    const TLorentzVector s = a + b;

                    if((*i)->VetoEnergy()==0 && (*j)->VetoEnergy()==0)
                        ggIM->Fill(s.M());

                    ttIM->Fill(s.M());
                }
                ++j;
            }
        }
        ++i;
    }
}


void CandidatesAnalysis::Finish()
{

}


void CandidatesAnalysis::ShowResult()
{
    canvas("CandidatesAnalysis")
            << ggIM << energy << theta << phi
            << nCandidatesEvent << CandMultiplicities
            << drawoption("colz")
            << cbdEE << cbtof << tapsdEE << tapstof << detectors
            << endc;
}

AUTO_REGISTER_PHYSICS(CandidatesAnalysis)
