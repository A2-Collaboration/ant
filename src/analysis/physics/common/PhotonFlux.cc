#include "PhotonFlux.h"

#include "slowcontrol/SlowControlVariables.h"
#include "base/Logger.h"
#include "expconfig/ExpConfig.h"
#include "base/PlotExt.h"

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

    TaggEff      = HistFac.makeTH1D("Tagg-Eff",              label_tagger, label_counts, bins_tagger,          "taggeff",       true);
    lifetime     = HistFac.makeTH1D("lifetime", "life time [%]",           label_counts, BinSettings(250,0,1), "lifetime",       true);
    lifetime_g   = HistFac.makeGraph("lifetime(t) [s]","lt_graph");

    ScalerCounts = HistFac.makeTH1D("scaler counts",         label_tagger, label_counts, bins_tagger,          "scalerCounts",  true);

    Flux         = HistFac.makeTH1D("Photonflux per channel",               label_tagger, label_counts,        bins_tagger, "flux",               true);
    FluxLTcor    = HistFac.makeTH1D("Photonflux per channel * lifetime",    label_tagger, label_counts,        bins_tagger, "fluxltcor",          true);

    IntLumi      = HistFac.makeTH1D("Integrated Luminosity", label_tagger, "Int. L [1/#mu b]",  bins_tagger, "intlumi",       true);
    IntLumiLTcor = HistFac.makeTH1D("Integrated Luminosity * lifetime", label_tagger, "Int. L [1/#mu b]",  bins_tagger, "intlumicor",       true);

    Lumi         = HistFac.makeTH1D("Luminosity",            label_tagger, "L [1/(#mu b * s)]", bins_tagger, "lumi",          true);

    info         = HistFac.makeTH1D("info", "","",BinSettings(2,0,0),"info");


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
    time += slowcontrol::Variables::Clocks->GetExpClock() / 1E6;
    const auto lt       = slowcontrol::Variables::Trigger->GetExpLivetime();
    GraphExt::FillGraph(lifetime_g,time,lt);

    for ( auto ch = 0u ; ch < nchannels ; ++ch)
    {
        const auto counts   = slowcontrol::Variables::TaggerScalers->GetCounts().at(ch);

        ScalerCounts->Fill(ch,counts);
        Flux->Fill(ch, counts);
        FluxLTcor->Fill(ch, counts * lt);
        lifetime->Fill(lt);
        IntLumi->Fill(ch, counts);
        IntLumiLTcor->Fill(ch, counts * lt);
        Lumi->Fill(ch,counts);
    }
}




void PhotonFlux::Finish()
{
    LOG(INFO) << "Total time recorded: " << time << " s";
    for ( auto ch = 0u ; ch < nchannels ; ++ch)
    {
        const auto taggEff  = tagger->GetTaggEff(ch);
        const auto bin = TaggEff->Fill(ch,taggEff.Value);
        TaggEff->SetBinError(bin, taggEff.Error);
    }

    Flux->Multiply(TaggEff);
    FluxLTcor->Multiply(TaggEff);

    IntLumi->Multiply(TaggEff);
    IntLumi->Scale(targetDensity);
    IntLumiLTcor->Multiply(TaggEff);
    IntLumiLTcor->Scale(targetDensity);

    Lumi->Multiply(TaggEff);
    Lumi->Scale(targetDensity);
    Lumi->Scale(1./time);
    Lumi->SetBit(TH1::kIsAverage);


    double intErr = 0;
    int    bin    = 0;
    bin = info->Fill("time [s]",time);
    info->SetBinError(bin,0);

    bin = info->Fill("total L^{int} [1/#mub]",IntLumi->IntegralAndError(1,nchannels,intErr));
    info->SetBinError(bin,intErr);

    LOG(INFO) << "Seen scaler-blocks: " << seenScalerBlocks;
}

void PhotonFlux::ShowResult()
{
    canvas("overview")
            << ScalerCounts
            << info
            << IntLumi
            << Lumi
            << endc;
}


AUTO_REGISTER_PHYSICS(PhotonFlux)
