#include "physics/common/DataOverview.h"
#include "plot/root_draw.h"
#include "expconfig/ExpConfig.h"
#include "root-addons/cbtaps_display/TH2CB.h"
#include "utils/combinatorics.h"
#include "base/Logger.h"

#include "TF1.h"

#include <cassert>

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace ant::analysis::data;

TaggerOverview::TaggerOverview(const string &name, PhysOptPtr opts):
    DataOverviewBase(name, opts)
{

    const auto tagger = ExpConfig::Setup::GetDetector<TaggerDetector_t>();
    const auto Emin   = tagger->GetPhotonEnergy(0);
    const auto Emax   = tagger->GetPhotonEnergy(tagger->GetNChannels()-1);
    const auto width  = (Emax - Emin)/tagger->GetNChannels();

    const BinSettings bins_hits(50);
    const BinSettings bins_energy(tagger->GetNChannels(), Emin-width/2, Emax+width/2);
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

void TaggerOverview::ProcessEvent(const Event &event)
{
    const auto taggerhits = (mode == Mode::Reconstructed) ? event.Reconstructed.TaggerHits : event.MCTrue.TaggerHits;

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


DataOverviewBase::DataOverviewBase(const string &name, PhysOptPtr opts):
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

const Event::Data &DataOverviewBase::GetBranch(const Event &event) const
{
   return (mode == Mode::Reconstructed) ? event.Reconstructed : event.MCTrue;
}


TriggerOverview::TriggerOverview(const string &name, PhysOptPtr opts):
    DataOverviewBase(name, opts)
{
    const BinSettings bins_errors(30);
    const BinSettings bins_multiplicity(10);
    const BinSettings bins_energy(1600);

    CBESum       = HistFac.makeTH1D("CB Energy Sum",  "CB Energy Sum [MeV]", "", bins_energy,      "CBESum");
    Multiplicity = HistFac.makeTH1D("Multiplicity",   "# Hits",              "", bins_multiplicity,"Multiplicity");
    nErrorsEvent = HistFac.makeTH1D("Errors / Event", "# errors",            "", bins_errors,      "nErrrorsEvent");


    triggereff = HistFac.make<TH2CB>("triggereff","Trigger Efficiency");

    auto cb_detector = ExpConfig::Setup::GetDetector(Detector_t::Type_t::CB);
    for(unsigned ch=0;ch<cb_detector->GetNChannels();ch++) {
        if(cb_detector->IsIgnored(ch)) {
            triggereff->CreateMarker(ch);
        }
    }

    triggereff_thetaphi = HistFac.makeTH2D("Trigger Efficiency Theta/Phi",
                                                 "#theta / #circ",
                                                 "#phi / #circ",
                                                 BinSettings(100,0,180),
                                                 BinSettings(150,-180,180),
                                                 "triggereff_thetaphi");
}

TriggerOverview::~TriggerOverview()
{}

void TriggerOverview::ProcessEvent(const Event &event)
{
    const auto& trigger = GetBranch(event).Trigger;

    CBESum->Fill(trigger.CBEnergySum);
    Multiplicity->Fill(trigger.ClusterMultiplicity);
    nErrorsEvent->Fill(trigger.Errors.size());


    // always run on Reconstructed
    const auto& cands = event.Reconstructed.Candidates;

    for( auto comb = analysis::utils::makeCombination(cands,2); !comb.Done(); ++comb ) {
        const CandidatePtr& p1 = comb.at(0);
        const CandidatePtr& p2 = comb.at(1);

        if(p1->VetoEnergy<0.25 && p2->VetoEnergy<0.25
           && (p1->Detector & Detector_t::Type_t::CB)
           && (p2->Detector & Detector_t::Type_t::CB)) {
            const Particle a(ParticleTypeDatabase::Photon,comb.at(0));
            const Particle b(ParticleTypeDatabase::Photon,comb.at(1));
            const TLorentzVector gg = a + b;

            auto cl1 = p1->FindCaloCluster();
            auto cl2 = p2->FindCaloCluster();
            const auto im = gg.M();
            if(im > 125 && im < 145) {
                triggereff->FillElement(cl1->CentralElement, 1.0);
                triggereff->FillElement(cl2->CentralElement, 1.0);
                triggereff_thetaphi->Fill(cl1->Position.Theta()*TMath::RadToDeg(),
                                          cl1->Position.Phi()*TMath::RadToDeg());
                triggereff_thetaphi->Fill(cl2->Position.Theta()*TMath::RadToDeg(),
                                          cl2->Position.Phi()*TMath::RadToDeg());
            }
        }
    }
}

void TriggerOverview::Finish()
{
    auto& h = triggereff_thetaphi;
    auto& h_proj = triggereff_proj;
    h_proj = h->ProjectionX("triggereff_proj");
    h_proj->Fit("pol1","","",30,155);
    TF1* func = h_proj->GetFunction("pol1");

    for(Int_t binx=0;binx<h->GetNbinsX()+1;binx++) {
        const auto theta = h->GetXaxis()->GetBinCenter(binx);
        const double corrfactor = func->Eval(theta) / func->Eval(90);
        for(Int_t biny=0;biny<h->GetNbinsY()+1;biny++) {
            if(corrfactor<0.1) {
                h->SetBinContent(binx, biny, 0);
                continue;
            }
            auto val = h->GetBinContent(binx, biny);
            h->SetBinContent(binx, biny, val/corrfactor);
        }
    }

    auto cb = ExpConfig::Setup::GetDetector(Detector_t::Type_t::CB);

    for(unsigned ch=0;ch<cb->GetNChannels();ch++) {
        if(cb->IsIgnored(ch))
            continue;
        const auto pos = cb->GetPosition(ch);
        const auto theta = pos.Theta()*TMath::RadToDeg();
        const double corrfactor = func->Eval(theta) / func->Eval(90);
        if(corrfactor<0.1) {
            triggereff->SetElement(ch, 0);
            continue;
        }
        auto val = triggereff->GetElement(ch);
        triggereff->SetElement(ch, val/corrfactor);
    }

    delete h_proj;
    h_proj = h->ProjectionX("triggereff_proj");
}

void TriggerOverview::ShowResult()
{
    canvas(this->GetName()+" "+GetMode())
            << CBESum
            << Multiplicity
            << nErrorsEvent
            << drawoption("colz") << triggereff
            << drawoption("colz") << triggereff_thetaphi
            << triggereff_proj
            << endc;
}

TargetOverview::TargetOverview(const string& name, PhysOptPtr opts):
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

void TargetOverview::ProcessEvent(const Event& event)
{
    const auto& target = GetBranch(event).Target;
    VertexXY->Fill(target.Vertex.X(), target.Vertex.Y());
    VertexZ->Fill(target.Vertex.Z());
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

ParticleOverview::ParticleOverview(const string &name, PhysOptPtr opts):
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

void ParticleOverview::ProcessEvent(const Event &event)
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

