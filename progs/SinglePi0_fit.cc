#include "base/Logger.h"

#include "tclap/CmdLine.h"
#include "base/interval.h"
#include "base/WrapTFile.h"
#include "base/std_ext/string.h"
#include "base/std_ext/system.h"
#include "base/std_ext/memory.h"
#include "base/std_ext/math.h"
#include "base/ParticleType.h"
#include "base/TH_ext.h"

#include "analysis/plot/RootDraw.h"
#include "root-addons/analysis_codes/Math.h"

#include "analysis/plot/HistogramFactory.h"

#include "TH1D.h"
#include "TH3D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPaveText.h"

#include "TSystem.h"
#include "TRint.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooArgusBG.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooHistPdf.h"
#include "RooHist.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "RooChebychev.h"
#include "RooConstVar.h"
#include "RooDerivative.h"
#include "RooFFTConvPdf.h"
#include "RooChi2Var.h"
#include "RooMinuit.h"
#include "RooFitResult.h"


using namespace ant;
using namespace std;
using namespace RooFit;

int main(int argc, char** argv) {
    SetupLogger();

    TCLAP::CmdLine cmd("SinglePi0_fit", ' ', "0.1");

    auto cmd_verbose       = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"int");

    auto cmd_data          = cmd.add<TCLAP::ValueArg<string>>("","data","Data input",true,"","rootfile");

    auto cmd_lumi          = cmd.add<TCLAP::ValueArg<string>>("","lumi","path to luminosity-class output",true,"","rootfile");

    auto cmd_eff           = cmd.add<TCLAP::ValueArg<string>>("","eff","MC signal input",true,"","rootfile");

    auto cmd_histpath      = cmd.add<TCLAP::ValueArg<string>>("","histpath","Path for hists (determines cutstr)",false,
                                                              "singlePi0_Plot/dicardedEk<20/EMB_prob>0.1/AllPhotonsInCB/NoTouchesHole/Pi0PIDVeto==0","path");
    auto cmd_histname      = cmd.add<TCLAP::ValueArg<string>>("","histname","Name of hist",false,"finalPlot","name");

    auto cmd_histluminame  = cmd.add<TCLAP::ValueArg<string>>("","histlumi","Name of hist",false,"intlumi","name");

    auto cmd_histreconame  = cmd.add<TCLAP::ValueArg<string>>("","histreco","Name of hist",false,"effrecon_pi0","name");
    auto cmd_histseenname  = cmd.add<TCLAP::ValueArg<string>>("","histseen","Name of hist",false,"seenMCcosTheta","name");



    auto cmd_batchmode = cmd.add<TCLAP::MultiSwitchArg>("b","batch","Run in batch mode (no ROOT shell afterwards)",false);
    auto cmd_output    = cmd.add<TCLAP::ValueArg<string>>("o","output","Output file",false,"","filename");


    cmd.parse(argc, argv);
    if(cmd_verbose->isSet()) {
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());
    }


    TH3D* h_data = nullptr;
    WrapTFileInput input_data(cmd_data->getValue()); // keep it open
    {
        const string histpath = cmd_histpath->getValue()+"/h/data/"+cmd_histname->getValue();
        if(!input_data.GetObject(histpath, h_data)) {
            LOG(ERROR) << "Cannot find " << histpath;
            return EXIT_FAILURE;
        }
    }

    TH1D* h_lumi= nullptr;
    WrapTFileInput input_lumi(cmd_lumi->getValue()); // keep it open
    {
        const string histpath = "PhotonFlux/"+cmd_histluminame->getValue();
        if(!input_lumi.GetObject(histpath, h_lumi)) {
            LOG(ERROR) << "Cannot find " << histpath;
            return EXIT_FAILURE;
        }
    }


    TH2D* h_eff= nullptr;
    WrapTFileInput input_eff(cmd_eff->getValue()); // keep it open
    {
        const string histpath_rec  = cmd_histpath->getValue() + "/h/Sig/" + cmd_histreconame->getValue();
        const string histpath_seen = "singlePi0_Plot/" + cmd_histseenname->getValue();
        TH2D* hseen = nullptr;

        if(!input_eff.GetObject(histpath_rec, h_eff)) {
            LOG(ERROR) << "Cannot find " << histpath_rec;
            return EXIT_FAILURE;
        }
        if(!input_eff.GetObject(histpath_seen, hseen)) {
            LOG(ERROR) << "Cannot find " << histpath_seen;
            return EXIT_FAILURE;
        }

        h_eff->Divide(hseen);
    }

    const auto nChannels       = h_data->GetNbinsY();
    const auto cosThetaBinning = TH_ext::getBins(h_data->GetXaxis());
//    const auto DeltaOmega      = cosThetaBinning.Length() * 2 * M_PI / cosThetaBinning.Bins();
//    const auto im2gBinning     = TH_ext::getBins(h_data->GetZaxis());


    // create TRint as RooFit internally creates functions/histograms, sigh...
    argc=0; // prevent TRint to parse any cmdline
    TRint app("OmegaEtaG_fit",&argc,argv,nullptr,0,true);

    unique_ptr<WrapTFileOutput> masterFile;
    if(cmd_output->isSet()) {
        // cd into masterFile upon creation
        masterFile = std_ext::make_unique<WrapTFileOutput>(cmd_output->getValue(), true);
    }

    analysis::HistogramFactory histfac("singlePi0_fits");
    vector<TH1D*> histChannels(nChannels);
    for (auto ch = 0 ; ch < nChannels ; ++ch)
    {

        histChannels.at(ch) = histfac.makeTH1D("Differential cross section",
                                              "cos(#theta_{coms})","diff. cross-section [#mu b / sr]",
                                              cosThetaBinning,std_ext::formatter() << "ch" << ch,true);
    }
    /*
    const interval<double> fitrange = TH_ext::getBins(h_data->GetXaxis());
    LOG(INFO) << "Fit Range: " << fitrange;
    const auto signalregion = interval<double>::CenterWidth(ParticleTypeDatabase::Omega.Mass(), 100.0);
    LOG(INFO) << "Signal Region: " << signalregion;


    // define observable and ranges
    RooRealVar var_IM("IM","IM", fitrange.Start(), fitrange.Stop(), "MeV");
    var_IM.setBins(10000);
    var_IM.setRange("full",var_IM.getMin(),var_IM.getMax());

    // load data to be fitted
    RooDataHist h_roo_data("h_roo_data","dataset",var_IM,h_data);

    // build shifted mc lineshape
    RooRealVar var_IM_shift("var_IM_shift", "shift in IM", 0.0, -10.0, 10.0);
    RooProduct var_IM_shift_invert("var_IM_shift_invert","shifted IM",RooArgSet(var_IM_shift, RooConst(-1.0)));
    RooAddition var_IM_shifted("var_IM_shifted","shifted IM",RooArgSet(var_IM,var_IM_shift_invert));
    RooDataHist h_roo_mc("h_roo_mc","MC lineshape", var_IM, h_mc);
    RooHistPdf pdf_mc_lineshape("pdf_mc_lineshape","MC lineshape as PDF", var_IM_shifted, var_IM, h_roo_mc, 2);

    // build gaussian
    RooRealVar  var_gauss_sigma("gauss_sigma","width of gaussian", 10.0, 0.0, 100.0);
    RooGaussian pdf_gaussian("pdf_gaussian","Gaussian smearing", var_IM, RooConst(0.0), var_gauss_sigma);

    // build signal as convolution, note that the gaussian must be the second PDF (see documentation)
    RooFFTConvPdf pdf_signal("pdf_signal","MC_lineshape (X) gauss",var_IM, pdf_mc_lineshape, pdf_gaussian) ;

    const int polOrder = 6;
    std::vector<std::unique_ptr<RooRealVar>> bkg_params; // RooRealVar cannot be copied, so create them on heap
    RooArgSet roo_bkg_params;
    for(int p=0;p<polOrder;p++) {
        bkg_params.emplace_back(std_ext::make_unique<RooRealVar>((
                                    "p_"+to_string(p)).c_str(), ("Bkg Par "+to_string(p)).c_str(), 0.0, -10.0, 10.0)
                                );
        roo_bkg_params.add(*bkg_params.back());
    }
    RooChebychev pdf_background("chebychev","Polynomial background",var_IM,roo_bkg_params);
    var_IM.setRange("bkg_l", var_IM.getMin(), signalregion.Start());
    var_IM.setRange("bkg_r", signalregion.Stop(), var_IM.getMax());
    pdf_background.fitTo(h_roo_data, Range("bkg_l,bkg_r")); // using Range(..., ...) does not work here (bug in RooFit, sigh)

    // build sum
    RooRealVar nsig("nsig","#signal events",     6E+6, 0, 1E+7);
    RooRealVar nbkg("nbkg","#background events", 9E+6, 0, 1E+9);
    RooAddPdf pdf_sum("pdf_sum","total sum",RooArgList(pdf_signal,pdf_background),RooArgList(nsig,nbkg));

    // do the actual maximum likelihood fit
    auto fr = pdf_sum.fitTo(h_roo_data, Extended(), SumW2Error(kTRUE), Range("full"), Save());
    const auto numParams = fr->floatParsFinal().getSize();

    // draw output, won't be shown in batch mode
    RooPlot* frame = var_IM.frame();
    h_roo_data.plotOn(frame);
    pdf_sum.plotOn(frame, LineColor(kRed));
    RooHist* hresid = frame->residHist();
    const auto chi2ndf = frame->chiSquare(numParams);
    pdf_sum.plotOn(frame, Components(pdf_background), LineColor(kBlue));
    pdf_sum.plotOn(frame, Components(pdf_signal), LineColor(kGreen));
    frame->Draw();

    new TCanvas();
    hresid->Draw();
    LOG(INFO) << "NSig: " << nsig.getValV();

    LOG(INFO) << "numParams=" << numParams << " chi2ndf=" << chi2ndf;

    fr->Print();
    */

    if(!cmd_batchmode->isSet()) {
        if(!std_ext::system::isInteractive()) {
            LOG(INFO) << "No TTY attached. Not starting ROOT shell.";
        }
        else {
            if(masterFile)
                LOG(INFO) << "Close ROOT properly to write data to disk.";

            app.Run(kTRUE); // really important to return...
            if(masterFile)
                LOG(INFO) << "Writing output file...";
            masterFile = nullptr;   // and to destroy the master WrapTFile before TRint is destroyed
        }
    }

    return 0;
}
