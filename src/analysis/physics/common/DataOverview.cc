#include "physics/common/DataOverview.h"
#include "plot/root_draw.h"
#include "expconfig/ExpConfig.h"
#include "root-addons/cbtaps_display/TH2CB.h"
#include "utils/combinatorics.h"
#include "base/Logger.h"

#include "TF1.h"
#include "TLorentzVector.h"

#include <cassert>

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

TaggerOverview::TaggerOverview(const string &name, OptionsPtr opts):
    DataOverviewBase(name, opts)
{

    const auto tagger = ExpConfig::Setup::GetDetector<TaggerDetector_t>();

    // it's not so trivial to generate a proper BinSettings avoiding binning artifacts
    // the following at least works nicely for the EPT
    vector<double> photon_energies(tagger->GetNChannels());
    for(unsigned ch=0;ch<tagger->GetNChannels();ch++)
        photon_energies[ch] = tagger->GetPhotonEnergy(ch);
    std::sort(photon_energies.begin(), photon_energies.end());

    const BinSettings bins_hits(50);
    const BinSettings bins_energy(BinSettings::Make(photon_energies));
    const BinSettings bins_channels(tagger->GetNChannels());
    const BinSettings bins_times(1000,-50,50);

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

TaggerOverview::~TaggerOverview()
{
}

void TaggerOverview::ProcessEvent(const TEvent& event, manager_t&)
{
    const auto taggerhits = (mode == Mode::Reconstructed) ? event.Reconstructed().TaggerHits : event.MCTrue().TaggerHits;

    for(auto hit=taggerhits.cbegin(); hit!=taggerhits.cend(); ++hit) {
        Channels->Fill(hit->Channel);
        Energies->Fill(hit->PhotonEnergy);
        Times->Fill(hit->Time);

        for(auto hit2 = next(hit); hit2!=taggerhits.cend(); ++hit2) {
            channel_correlation->Fill(hit->Channel, hit2->Channel);
            channel_correlation->Fill(hit2->Channel, hit->Channel);
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


DataOverviewBase::DataOverviewBase(const string &name, OptionsPtr opts):
    Physics(name, opts)
{
    if(opts->Get<string>("Mode") == "Reconstructed")
        mode = Mode::Reconstructed;
    else if(opts->Get<string>("Mode") == "MCTrue")
        mode = Mode::MCTrue;

    HistFac.SetTitlePrefix(GetMode());
}

DataOverviewBase::~DataOverviewBase()
{}

string DataOverviewBase::GetMode() const
{
    if(mode == Mode::Reconstructed) {
        return "Reconstructed";
    } else
        return "MCTrue";
}

const TEventData& DataOverviewBase::GetBranch(const TEvent& event) const
{
   return (mode == Mode::Reconstructed) ? event.Reconstructed() : event.MCTrue();
}


TriggerOverview::TriggerOverview(const string &name, OptionsPtr opts):
    DataOverviewBase(name, opts)
{
    const BinSettings bins_errors(30);
    const BinSettings bins_multiplicity(10);
    const BinSettings bins_energy(1600);

    CBESum       = HistFac.makeTH1D("CB Energy Sum",  "CB Energy Sum [MeV]", "", bins_energy,      "CBESum");
    Multiplicity = HistFac.makeTH1D("Multiplicity",   "# Hits",              "", bins_multiplicity,"Multiplicity");
    nErrorsEvent = HistFac.makeTH1D("Errors / Event", "# errors",            "", bins_errors,      "nErrorsEvent");
    CBTiming     = HistFac.makeTH1D("CB Timing", "CB Timing [ns]", "", BinSettings(300,-20,20), "CBTiming");


    auto cb_detector = ExpConfig::Setup::GetDetector(Detector_t::Type_t::CB);

    CBESum_perCh = HistFac.makeTH2D("CBEsum vs. CB channel",
                                    "CBEsum / MeV",
                                    "Channel",
                                    bins_energy,
                                    BinSettings(cb_detector->GetNChannels()),
                                    "CBESum_perCh");
    E_perCh = HistFac.makeTH2D("E vs. CB channel",
                               "E / MeV",
                               "Channel",
                               bins_energy,
                               BinSettings(cb_detector->GetNChannels()),
                               "E_perCh");
    t = HistFac.makeTTree("trigger");

    tree.CreateBranches(t);

}

TriggerOverview::~TriggerOverview()
{}

void TriggerOverview::ProcessEvent(const TEvent& event, manager_t&)
{
    const auto& branch =  GetBranch(event);
    const auto& trigger = branch.Trigger;

    CBESum->Fill(trigger.CBEnergySum);
    tree.CBESum = trigger.CBEnergySum;

    Multiplicity->Fill(trigger.ClusterMultiplicity);
    nErrorsEvent->Fill(trigger.DAQErrors.size());
    CBTiming->Fill(trigger.CBTiming);

    for(const TCluster& cluster : branch.Clusters) {
        if(cluster.DetectorType == Detector_t::Type_t::CB) {
            for(const TClusterHit& hit : cluster.Hits) {
                CBESum_perCh->Fill(trigger.CBEnergySum, hit.Channel);
                E_perCh->Fill(hit.Energy, hit.Channel);
            }
        }
    }

    const auto& protons = event.MCTrue().Particles.Get(ParticleTypeDatabase::Proton);
    if(!protons.empty()) {
        tree.true_proton() = *protons.front();
    } else {
        tree.true_proton() = TLorentzVector();
    }

    for(const auto& th : branch.TaggerHits) {
        tree.TaggE = th.PhotonEnergy;
        tree.TaggT = th.Time;
        t->Fill();
    }
}

void TriggerOverview::Finish()
{

}

void TriggerOverview::ShowResult()
{
    canvas(this->GetName()+" "+GetMode())
            << CBESum
            << Multiplicity
            << nErrorsEvent
            << drawoption("colz") << CBESum_perCh
            << E_perCh
            << CBTiming
            << endc;
}

TargetOverview::TargetOverview(const string& name, OptionsPtr opts):
    DataOverviewBase(name, opts)
{
    VertexXY = HistFac.makeTH2D("Vertex XY","X/cm","Y/cm",
                                BinSettings(100,-5,5), BinSettings(100,-5,5),
                                "VertexXY");
    VertexZ = HistFac.makeTH1D("Vertex Z","Z/cm","#",
                               BinSettings(400,-10,10),
                               "VertexZ");
}

TargetOverview::~TargetOverview()
{}

void TargetOverview::ProcessEvent(const TEvent& event, manager_t&)
{
    const auto& target = GetBranch(event).Target;
    VertexXY->Fill(target.Vertex.x, target.Vertex.y);
    VertexZ->Fill(target.Vertex.z);
}

void TargetOverview::ShowResult()
{
    canvas(this->GetName()+" "+GetMode())
            << drawoption("colz") << VertexXY
            << VertexZ
            << endc;
}


void ParticleOverview::SetBinLabels(TH1D *hist, const ParticleTypeDatabase::TypeList_t& types)
{
    assert(int(types.size()) <= hist->GetNbinsX());
    int i=1;
    for(const auto& type : types) {
        hist->GetXaxis()->SetBinLabel(i++, type->PrintName().c_str());
    }
}

ParticleOverview::ParticleOverview(const string &name, OptionsPtr opts):
    DataOverviewBase(name, opts)
{
    const BinSettings bins_particles(15);

    nParticles    = HistFac.makeTH1D("Particles / Event", "# Particles", "", bins_particles, "nParticles");
    particleTypes = HistFac.makeTH1D("Particle Types",    "Type",        "", BinSettings(unsigned(ParticleTypeDatabase::DetectableTypes().size())), "ParticleTypes");
    SetBinLabels(particleTypes, ParticleTypeDatabase::DetectableTypes());

    for(const ParticleTypeDatabase::Type* type : ParticleTypeDatabase::DetectableTypes()) {
        nType[type] = HistFac.makeTH1D(type->PrintName() + " / Event", "# Particles", "", bins_particles, type->PrintName());
    }
}

ParticleOverview::~ParticleOverview()
{

}

void ParticleOverview::ProcessEvent(const TEvent& event, manager_t&)
{
    const auto& particles = GetBranch(event).Particles.GetAll();
    nParticles->Fill(particles.size());

    for(const auto& p : particles) {
        particleTypes->Fill(p->Type().PrintName().c_str(), 1.0);
    }

    for(const ParticleTypeDatabase::Type* type : ParticleTypeDatabase::DetectableTypes()) {
        nType[type]->Fill(GetBranch(event).Particles.Get(*type).size());
    }
}

void ParticleOverview::ShowResult()
{
    canvas c(GetName());

    c << nParticles
      << particleTypes;

    for(const auto& entry : nType) {
        c<< entry.second;
    }

    c  << endc;
}

AUTO_REGISTER_PHYSICS(ParticleOverview)
AUTO_REGISTER_PHYSICS(TaggerOverview)
AUTO_REGISTER_PHYSICS(TriggerOverview)
AUTO_REGISTER_PHYSICS(TargetOverview)

