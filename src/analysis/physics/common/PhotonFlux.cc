#include "PhotonFlux.h"

#include "slowcontrol/SlowControlVariables.h"
#include "base/Logger.h"
#include "expconfig/ExpConfig.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

PhotonFlux::PhotonFlux(const std::string& name, OptionsPtr opts) :
    Physics(name, opts),
    tagger(ExpConfig::Setup::GetDetector<TaggerDetector_t>()),
    nchannels(tagger->GetNChannels())
{
    const auto bins_tagger    = BinSettings(nchannels);
    const string label_tagger = "tagger channel";
    const string label_counts = "#";
    ScalerCounts = HistFac.makeTH1D("scaler counts", label_tagger, label_counts, bins_tagger, "scalerCounts", true);
    Flux         = HistFac.makeTH1D("Photon Flux",   label_tagger, label_counts, bins_tagger, "flux",         true);
    Lumi         = HistFac.makeTH1D("Luminosity",    label_tagger, label_counts, bins_tagger, "lumi",         true);
    TaggEff      = HistFac.makeTH1D("Tagg-Eff",      label_tagger, label_counts, bins_tagger, "taggeff",      true);

    slowcontrol::Variables::TaggerScalers->Request();
    slowcontrol::Variables::Clocks->Request();
    slowcontrol::Variables::Trigger->Request();
    slowcontrol::Variables::Beam->Request();
}


PhotonFlux::~PhotonFlux() {}

void PhotonFlux::ProcessEvent(const TEvent&, manager_t& )
{
    if(slowcontrol::Variables::TaggerScalers->HasChanged())
    {
        seenScalerBlocks++;
        processBlock();
    }
}


void PhotonFlux::processBlock()
{
    for ( auto ch = 0u ; ch < nchannels ; ++ch)
    {
        const auto counts   = slowcontrol::Variables::TaggerScalers->GetCounts().at(ch);
        ScalerCounts->Fill(ch,counts);
        Flux->Fill(ch, counts);
        Lumi->Fill(ch, counts);
    }
}




void PhotonFlux::Finish()
{
    for ( auto ch = 0u ; ch < nchannels ; ++ch)
    {
        const auto taggEff  = tagger->GetTaggEff(ch);
        const auto bin = TaggEff->Fill(ch,taggEff.Value);
        TaggEff->SetBinError(bin, taggEff.Error);
    }
    Flux->Divide(TaggEff);
    Lumi->Divide(TaggEff);
    Lumi->Scale(targetDensity);
    LOG(INFO) << "Seen scaler-blocks: " << seenScalerBlocks;
}

void PhotonFlux::ShowResult()
{
    canvas("overview")
            << ScalerCounts
            << Flux
            << Lumi
            << endc;
}


AUTO_REGISTER_PHYSICS(PhotonFlux)
