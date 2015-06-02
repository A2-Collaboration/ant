#include "analysis/DataOverview.h"
#include "plot/root_draw.h"

using namespace ant;
using namespace std;


analysis::DataOverview::OverviewSet::OverviewSet(SmartHistFactory &factory, const string &title)
{
    SmartHistFactory myfac(title, factory);

    TaggerChannel = myfac.makeHist<int>("Tagger Channel ("+title+")","Channel","",BinSettings(400),"taggerchannel");
    PhotonEnergy  = myfac.makeHist<double>("Tagger Spectrum (Photon Energy) ("+title+")","E_{Beam} [MeV]","",BinSettings(2000),"taggerspectrum");
    TaggedTime    = myfac.makeHist<double>("Tagged Time ("+title+")","t [ns]","",BinSettings(1000,-100,100),"taggedtime");

    nParticles    = myfac.makeHist<int>("Number of particles/event ("+title+")","Particles/Event","",BinSettings(16),"nParticles");
    ParticleTypes = myfac.makeHist<string>("Particle Types ("+title+")","Particle Type","#",BinSettings(1),"ParticleTypes");

    CBEnergySum   = myfac.makeHist<double>("CB Energy Sum ("+title+")","E_{CB} [MeV]","",BinSettings(1600),"CBEsum");


}

void analysis::DataOverview::OverviewSet::Fill(const Event::Data &dataset)
{
    for(auto& taggerhit : dataset.TaggerHits() ) {
        TaggerChannel.Fill(taggerhit->Channel());
        PhotonEnergy.Fill(taggerhit->PhotonEnergy());
        TaggedTime.Fill(taggerhit->Time());
    }

    nParticles.Fill(dataset.Particles().GetAll().size());
    for(auto& particle : dataset.Particles().GetAll()) {
        ParticleTypes.Fill(particle->Type().PrintName());
    }

    CBEnergySum.Fill(dataset.TriggerInfos().CBEenergySum());
}



analysis::DataOverview::DataOverview(const string &name):
    Physics(name),
    reconstructed(HistFac, "reconstructed"),
    mctrue(HistFac,"mctrue")
{
}

void ant::analysis::DataOverview::ProcessEvent(const ant::Event &event)
{
    reconstructed.Fill(event.Reconstructed());
    mctrue.Fill(event.MCTrue());
}


void ant::analysis::DataOverview::Finish()
{

}


void ant::analysis::DataOverview::ShowResult()
{
    canvas("DataOverview: Tagger")
            << reconstructed.TaggerChannel << mctrue.TaggerChannel
            << reconstructed.TaggedTime    << mctrue.TaggedTime
            << reconstructed.PhotonEnergy  << mctrue.PhotonEnergy
            << endc;

    canvas("DataOverview: Particle Stats")
            << reconstructed.nParticles    << mctrue.nParticles
            << reconstructed.ParticleTypes << mctrue.ParticleTypes
            << endc;

    canvas("DataOverview: TriggerInfo")
            << reconstructed.CBEnergySum   << mctrue.CBEnergySum
            << endc;
}


