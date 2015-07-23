#include "physics/common/TrackAnalysis.h"
#include "plot/root_draw.h"
#include "data/Particle.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include <cmath>

using namespace std;
using namespace ant;
using namespace ant::analysis;


TrackAnalysis::TrackAnalysis(const string &name):
    Physics(name)
{
    nTracksEvent = HistFac.makeTH1D("nTracks/Event","Tracks","",BinSettings(30),"nTracks");
    energy = HistFac.makeTH1D("Energy","E [MeV]","",BinSettings(1000),"energy");
    theta  = HistFac.makeTH1D("Theta","#theta [#circ]","",BinSettings(360,0,180),"theta");
    phi    = HistFac.makeTH1D("Phi","#phi [#circ]","",BinSettings(720,-180,180),"phi");
    ggIM         = HistFac.makeTH1D("2 Neutral Tracks IM","M [MeV]","",BinSettings(1000),"ggIM");
}

void TrackAnalysis::ProcessEvent(const ant::Event &event)
{
    const auto& tracks = event.Reconstructed().Tracks();

    nTracksEvent->Fill(tracks.size());

    TrackList::const_iterator i = tracks.begin();

    while(i!=tracks.end()) {
        energy->Fill((*i)->ClusterEnergy());
        theta->Fill((*i)->Theta()*TMath::RadToDeg());
        phi->Fill((*i)->Phi()*TMath::RadToDeg());

        if(std::isnan((*i)->VetoEnergy()) && (*i)->ClusterEnergy()>20.0) {
            const Particle a(ParticleTypeDatabase::Photon,*i);
            TrackList::const_iterator j = i;
            ++j;
            while(j!=tracks.end()) {
                if(std::isnan((*j)->VetoEnergy())&& (*i)->ClusterEnergy()>20.0) {
                    const Particle b(ParticleTypeDatabase::Photon,*j);
                    const TLorentzVector s = a + b;
                    ggIM->Fill(s.M());
                }
                ++j;
            }
        }
        ++i;
    }
}


void TrackAnalysis::Finish()
{

}


void TrackAnalysis::ShowResult()
{
    canvas("TrackAnalysis")
            << ggIM << energy << theta << phi
            << nTracksEvent
            << endc;
}


