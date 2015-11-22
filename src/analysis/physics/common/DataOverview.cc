#include "physics/common/DataOverview.h"
#include "plot/root_draw.h"
#include "expconfig/ExpConfig.h"


using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace ant::analysis::data;

DataOverview::OverviewSet::OverviewSet(SmartHistFactory &factory, const string &title)
{
    SmartHistFactory myfac(title, factory);

    TaggerChannel = myfac.makeHist<int>("Tagger Channel ("+title+")","Channel","",BinSettings(400),"taggerchannel");
    PhotonEnergy  = myfac.makeHist<double>("Tagger Spectrum (Photon Energy) ("+title+")","E_{Beam} [MeV]","",BinSettings(2000),"taggerspectrum");
    TaggedTime    = myfac.makeHist<double>("Tagged Time ("+title+")","t [ns]","",BinSettings(1000,-100,100),"taggedtime");

    nParticles    = myfac.makeHist<int>("Number of particles/event ("+title+")","Particles/Event","",BinSettings(16),"nParticles");
    ParticleTypes = myfac.makeHist<string>("Particle Types ("+title+")","Particle Type","#",BinSettings(1),"ParticleTypes");

    CBEnergySum   = myfac.makeHist<double>("CB Energy Sum ("+title+")","E_{CB} [MeV]","",BinSettings(1600),"CBEsum");


}

void DataOverview::OverviewSet::Fill(const Event::Data &dataset)
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



DataOverview::DataOverview(const std::string& name, PhysOptPtr opts):
    Physics(name, opts),
    reconstructed(HistFac, "reconstructed"),
    mctrue(HistFac,"mctrue")
{
}

void DataOverview::ProcessEvent(const Event &event)
{
    reconstructed.Fill(event.Reconstructed());
    mctrue.Fill(event.MCTrue());
}


void DataOverview::Finish()
{

}


void DataOverview::ShowResult()
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



string TaggerOverview::GetMode() const
{
    if(mode == Mode::Reconstructed) {
        return "Reconstructed";
    } else
        return "MCTrue";
}

TaggerOverview::TaggerOverview(const string &name, PhysOptPtr opts):
    Physics(name, opts)
{

    if(opts->GetOption("Mode") == "Reconstructed")
        mode = Mode::Reconstructed;
    else if(opts->GetOption("Mode") == "MCTrue")
        mode = Mode::MCTrue;

    HistFac.SetTitlePrefix(GetMode());

    const BinSettings bins_hits(50);
    const BinSettings bins_energy(200, 1400, 1600);
    const BinSettings bins_channels(47);
    const BinSettings bins_times(1000,-50,50);

//    const auto setup = ExpConfig::Setup::GetLastFound();
//    const auto tagger = setup ? setup->GetDetector<TaggerDetector_t>() : nullptr;

//    if(!tagger) {
//        throw std::runtime_error("No Tagger in Setup!");
//    }

//    const BinSettings bins_channels(tagger->GetNChannels());

//    const auto Emin   = tagger->GetPhotonEnergy(0);
//    const auto Emax   = tagger->GetPhotonEnergy(tagger->GetNChannels()-1);
//    const auto width  = (Emax - Emin)/tagger->GetNChannels();

//    const BinSettings bins_energy(tagger->GetNChannels(), Emin-width/2, Emax+width/2);

    nHitsEvent = HistFac.makeTH1D("Tagger Hits / Enent ", "# Hits/Event",   "",     bins_hits,    "TaggerHitsPerEvent");
    nHitsEvent->SetFillColor(kRed);

    Channels   = HistFac.makeTH1D("Tagger Channels hit ", "Channel Number", "Hits", bins_channels," TaggerChannels");
    Channels->SetFillColor(kYellow);

    Energies   = HistFac.makeTH1D("Tagged Photon Energies ", "E_{#gamma,tag} [MeV]", "", bins_energy, "TaggedEnergies");
    Energies->SetFillColor(kBlue);

    Times      = HistFac.makeTH1D("Tagger Hit Times ", "Time [ns]", "", bins_times, "TaggedTimes");
    Times->SetFillColor(kGreen);

    channel_correlation = HistFac.makeTH2D("Tagger Channel Correlation","Channel", "Channel", bins_channels, bins_channels, "ChannelCorreleation");
}

void TaggerOverview::ProcessEvent(const Event &event)
{
    const auto taggerhits = (mode == Mode::Reconstructed) ? event.Reconstructed().TaggerHits() : event.MCTrue().TaggerHits();

    for(auto hit=taggerhits.cbegin(); hit!=taggerhits.cend(); ++hit) {
        Channels->Fill((*hit)->Channel());
        Energies->Fill((*hit)->PhotonEnergy());
        Times->Fill((*hit)->Time());

        for(auto hit2 = next(hit); hit2!=taggerhits.cend(); ++hit2) {
            channel_correlation->Fill((*hit)->Channel(), (*hit2)->Channel());
            channel_correlation->Fill((*hit2)->Channel(), (*hit)->Channel());
        }
    }
    nHitsEvent->Fill(taggerhits.size());
}

void TaggerOverview::ShowResult()
{
    canvas(this->GetName()+" "+GetMode())
            << nHitsEvent
            << Channels
            << Energies
            << Times
            << drawoption("colz") << channel_correlation
            << endc;

}


AUTO_REGISTER_PHYSICS(DataOverview)
AUTO_REGISTER_PHYSICS(TaggerOverview)

