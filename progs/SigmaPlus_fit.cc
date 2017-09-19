#include "base/Logger.h"
#include "detail/SigmaPlus_tools.h"

#include "tclap/CmdLine.h"
#include "tclap/ValuesConstraintExtra.h"
#include "base/WrapTFile.h"
#include "analysis/utils/ValError.h"
#include "expconfig/ExpConfig.h"

#include "analysis/plot/HistogramFactory.h"
#include "analysis/utils/TaggerBins.h"

#include "base/ParticleType.h"



#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TGraphErrors.h"

#include "TSystem.h"
#include "TRint.h"


using namespace ant;
using namespace std;
using namespace RooFit;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::utils;
using namespace SIGMA;


int main(int argc, char** argv) {
    SetupLogger();

    TCLAP::CmdLine cmd("SinglePi0_fit", ' ', "0.1");

    auto cmd_verbose       = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"int");

    auto cmd_PlotFile      = cmd.add<TCLAP::ValueArg<string>>("","plotFile","Output from sigmaPlus_FinalPlot", true,"","rootfile");
    auto cmd_relpath_data  = cmd.add<TCLAP::ValueArg<string>>("","relPath", "Path to inside folders for bins",               false, "h/Data/","string");

    auto cmd_MCTrue        = cmd.add<TCLAP::ValueArg<string>>("","mcFile","Output from sigmaPlus_FinalPlot for MC True sample", true,"","rootfile");
    auto cmd_relpath_mc    = cmd.add<TCLAP::ValueArg<string>>("","relPathMC", "Path to inside folders for bins",               false, "h/Sig/","string");

    auto cmd_HistPath      = cmd.add<TCLAP::ValueArg<string>>("","histPath","Path to histgrams",               false, "sigmaPlus_FinalPlot/","string");
    auto cmd_HistName      = cmd.add<TCLAP::ValueArg<string>>("","histName","Data input from  singlePi0-Plot", false,"ppi0_2pi0","rootfile");

    auto cmd_nEgBins       = cmd.add<TCLAP::ValueArg<unsigned>>("","nEgBins","Number of bins in Egamma (only override if you know what you are doing)",false,10,"unsigned");

    auto cmd_lumi          = cmd.add<TCLAP::ValueArg<string>>("","lumi","path to luminosity-class output if seperate file from data input",true,"","rootfile");
    auto cmd_histluminame  = cmd.add<TCLAP::ValueArg<string>>("","histlumi","Name of hist",false,"intlumicor","name");
    auto cmd_histreconame  = cmd.add<TCLAP::ValueArg<string>>("","histreco","Name of hist",false,"recon_fit","name");
    auto cmd_seenprefix    = cmd.add<TCLAP::ValueArg<string>>("","seenprefix","Name of hist",false,"hist00","name");

    auto cmd_output        = cmd.add<TCLAP::ValueArg<string>>("o","output","Output file",false,"","filename");


    TCLAP::ValuesConstraintExtra<decltype(ExpConfig::Setup::GetNames())> allowedsetupnames(ExpConfig::Setup::GetNames());
    auto cmd_setup = cmd.add<TCLAP::ValueArg<string>>("s","setup","Setup to determine tagged photon energy bins",false,"Setup_2014_10_EPT_Prod",&allowedsetupnames);



    cmd.parse(argc, argv);
    if(cmd_verbose->isSet()) {
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());
    }
    //    auto loadHist = [](const WrapTFileInput& input, const string& histpath)
    //    {
    //        TH1* hist = nullptr;
    //        if (!input.GetObject(histpath,hist))
    //            throw runtime_error(std_ext::formatter() << "Cannot find " << histpath);
    //        return hist;
    //    };


    ExpConfig::Setup::SetByName(cmd_setup->getValue());

    argc=0;
    TRint app("SigmaPlus_fit",&argc,argv,nullptr,0,true);

    unique_ptr<WrapTFileOutput> outFile;
    if(cmd_output->isSet()) {
        // cd into masterFile upon creation
        outFile = std_ext::make_unique<WrapTFileOutput>(cmd_output->getValue(), true);
    }
    HistogramFactory histfac("SigmaPlus_fit");

    // get photonlux histogram
    const string fluxPath    = "PhotonFlux/";
    WrapTFileInput input_lumi(cmd_lumi->getValue());
    auto tagger = ExpConfig::Setup::GetDetector<TaggerDetector_t>();
    auto h_lumi = functions::getHist<TH1D>(input_lumi, fluxPath + cmd_histluminame->getValue());
    auto h_lumi_Eg = TaggerBins::ChannelToEg(h_lumi,tagger);


    LumiFitter_t fitter;
    const auto fitResult_Eg = fitter.DoFit(h_lumi_Eg);
    const auto fitResult = fitter.DoFit(h_lumi);
    auto fitfkt = LumiFitter_t::makeROOTfunction(fitResult,0,47);
    auto fitfkt_Eg = LumiFitter_t::makeROOTfunction(fitResult_Eg,1420,1580);

    //tagger-Energies
    const auto nEgammaBins = cmd_nEgBins->getValue();
//    const auto taggerBinRanges = utils::TaggerBins::MakeEgBins(tagger,nEgammaBins);
    const auto taggerBinRanges = utils::TaggerBins::EPTBinning();


    //     auto histEff = histfac.makeTH2D("eff",
    //                                     taggerLabel,cosThetaLabel,
    //                                     taggBins,cosThetaBins,"eff",true);
    //     histEff->Add(h_rec);
    //     histEff->Divide(h_);


    LOG(INFO) << "Loading data hists..." ;
    WrapTFileInput input_plotter(cmd_PlotFile->getValue());
    auto dalitzHists = functions::getHists(input_plotter,cmd_HistPath->getValue(), cmd_relpath_data->getValue(), cmd_HistName->getValue());
    LOG(INFO) << "Loading mc hists..." ;
    WrapTFileInput input_plotter_mc(cmd_MCTrue->getValue());
    auto dalitzHists_MC = functions::getHists(input_plotter_mc, cmd_HistPath->getValue(), cmd_relpath_mc->getValue(), cmd_HistName->getValue());
    const auto seenMCs = functions::getSeenMCsums(input_plotter_mc, cmd_HistPath->getValue(), cmd_seenprefix->getValue());

    auto dataGraph = histfac.makeGraphErrors("counts data","nData");
    auto effGraph = histfac.makeGraphErrors("counts mc","nMC");
    auto finalGraph = histfac.makeGraphErrors("cross section","finalPlot");
    auto lumiGraph = histfac.makeGraphErrors("int. luminosity","lumiGraph");

    const auto branchingRatio =
            BR::is_K0S
            * BR::Sigma_Pi0p * BR::Pi0_gg
            * BR::K0S_Pi0Pi0 * BR::Pi0_gg * BR::Pi0_gg;

    for (auto i = 0u ; i < TaggerBins::EPTBinning().size() ; ++i)
    {
        const auto e            = taggerBinRanges.at(i);

        const auto result       = tools::fitHist(dalitzHists.at(i));
        const auto resultmc     = tools::fitHist(dalitzHists_MC.at(i));
        const auto eff  = resultmc / seenMCs.at(i);
        const auto resultEffCor = result / eff;
        const auto integralLumi = ValError::Statistical(fitfkt->Integral(e.Start(),e.Stop()));
        const auto sigma_full  = resultEffCor  / integralLumi / branchingRatio;

        if (!isnan(result.v) && !isnan(resultmc.v) )
        {
            tools::FillGraphErrors(dataGraph, e.Center(),result.v, 0 ,result.e);
            tools::FillGraphErrors(effGraph, e.Center(),eff.v, 0 ,eff.e);
            tools::FillGraphErrors(finalGraph, e.Center(),sigma_full.v, 0 ,sigma_full.e);
            tools::FillGraphErrors(lumiGraph, e.Center(),integralLumi.v, 0 ,integralLumi.e);
        }
    }


    auto cdata = new TCanvas("result data","result data",600,600);
    cdata->Divide(2,2);
    cdata->cd(1);
    dataGraph->Draw("AP");

    cdata->cd(2);
    effGraph->Draw("AP");

    cdata->cd(3);
    finalGraph->Draw("AP");

    cdata->cd(4);
    h_lumi_Eg->Draw("");
    lumiGraph->Draw("AP same");
    fitfkt_Eg->Draw("same");

    {
        auto setLabels = [](TGraphErrors* g, const string& y)
        {
            g->GetXaxis()->SetTitle("E_{#gamma} [MeV]");
            g->GetYaxis()->SetTitle(y.c_str());
        };
        setLabels(dataGraph,"#");
        setLabels(effGraph,"eff [%]");
        setLabels(finalGraph,"#sigma [#mub]");
    }

    if(outFile)
        LOG(INFO) << "Close ROOT properly to write data to disk.";

    app.Run(kTRUE); // really important to return...
    if(outFile)
        LOG(INFO) << "Writing output file...";
    outFile = nullptr;

    return 0;
}
