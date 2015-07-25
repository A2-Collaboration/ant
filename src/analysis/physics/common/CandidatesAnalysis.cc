#include "physics/common/CandidatesAnalysis.h"
#include "plot/root_draw.h"
#include "data/Particle.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include <cmath>

using namespace std;
using namespace ant;
using namespace ant::analysis;


CandidatesAnalysis::CandidatesAnalysis(const string &name):
    Physics(name)
{
    nCandidatesEvent = HistFac.makeTH1D("Candidates/Event","Candidates","",BinSettings(30),"nCand");
    energy = HistFac.makeTH1D("Energy","E [MeV]","",BinSettings(1000),"energy");
    theta  = HistFac.makeTH1D("Theta","#theta [#circ]","",BinSettings(360,0,180),"theta");
    phi    = HistFac.makeTH1D("Phi","#phi [#circ]","",BinSettings(720,-180,180),"phi");
    ggIM         = HistFac.makeTH1D("2 Neutral Candidates IM","M [MeV]","",BinSettings(1000),"ggIM");
    ttIM         = HistFac.makeTH1D("2 Candidates IM","M [MeV]","",BinSettings(1000),"ttIM");
}

void CandidatesAnalysis::ProcessEvent(const ant::Event &event)
{
    const auto& candidates = event.Reconstructed().Candidates();

    nCandidatesEvent->Fill(candidates.size());

    CandidateList::const_iterator i = candidates.begin();

    while(i!=candidates.end()) {
        energy->Fill((*i)->ClusterEnergy());
        theta->Fill((*i)->Theta()*TMath::RadToDeg());
        phi->Fill((*i)->Phi()*TMath::RadToDeg());

        if((*i)->ClusterEnergy()>20.0) {
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
            << nCandidatesEvent
            << endc;
}

AUTO_REGISTER_PHYSICS(CandidatesAnalysis, "CandidatesAnalysis")
