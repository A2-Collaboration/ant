#include "base/Logger.h"

#include "tclap/CmdLine.h"
#include "base/interval.h"
#include "base/WrapTFile.h"
#include "base/std_ext/string.h"
#include "base/std_ext/system.h"
#include "base/std_ext/memory.h"
#include "base/ParticleType.h"

#include "analysis/plot/RootDraw.h"
#include "root-addons/analysis_codes/Math.h"

#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TLatex.h"

#include "TSystem.h"
#include "TRint.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooArgusBG.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooHistPdf.h"
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
#include "RooHist.h"
#include "RooPlotable.h"
#include "RooTrace.h"
#include "RooCBShape.h"

using namespace ant;
using namespace std;
using namespace RooFit;

int main(int argc, char** argv) {
    RooTrace::active(kTRUE) ;
    RooTrace::verbose(kTRUE) ;

    SetupLogger();

    TCLAP::CmdLine cmd("EtapOmegaG_fit", ' ', "0.1");
    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"int");
    auto cmd_data = cmd.add<TCLAP::ValueArg<string>>("","data","Data input",true,"","rootfile");
    auto cmd_mc = cmd.add<TCLAP::ValueArg<string>>("","mc","MC signal/reference input",true,"","rootfile");
    auto cmd_histpath = cmd.add<TCLAP::ValueArg<string>>("","histpath","Path for hists (determines cutstr)",false,"EtapOmegaG_plot_Ref/DiscardedEk=0/KinFitProb>0.02","path");
    auto cmd_histname = cmd.add<TCLAP::ValueArg<string>>("","histname","Name of hist",false,"h_IM_2g","name");
    auto cmd_batchmode = cmd.add<TCLAP::MultiSwitchArg>("b","batch","Run in batch mode (no ROOT shell afterwards)",false);
    auto cmd_output = cmd.add<TCLAP::ValueArg<string>>("o","output","Output file",false,"","filename");


    cmd.parse(argc, argv);
    if(cmd_verbose->isSet()) {
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());
    }


    TH1D* h_data = nullptr;
    WrapTFileInput input_data(cmd_data->getValue()); // keep it open
    {
        const string histpath = cmd_histpath->getValue()+"/h/Data/"+cmd_histname->getValue();
        if(!input_data.GetObject(histpath, h_data)) {
            LOG(ERROR) << "Cannot find " << histpath;
            return EXIT_FAILURE;
        }
    }

    TH1D* h_mc = nullptr;
    WrapTFileInput input_mc(cmd_mc->getValue()); // keep it open
    {
        const string histpath = cmd_histpath->getValue()+"/h/Sum_MC/"+cmd_histname->getValue();
        if(!input_mc.GetObject(histpath, h_mc)) {
            LOG(ERROR) << "Cannot find " << histpath;
            return EXIT_FAILURE;
        }
    }

    // create TRint as RooFit internally creates functions/histograms, sigh...
    argc=0; // prevent TRint to parse any cmdline
    TRint app("EtapOmegaG_plot",&argc,argv,nullptr,0,true);

    unique_ptr<WrapTFileOutput> masterFile;
    if(cmd_output->isSet()) {
        // cd into masterFile upon creation
        masterFile = std_ext::make_unique<WrapTFileOutput>(cmd_output->getValue(), true);
    }

    const interval<double> signal_region{920, 990};

    // define observable and ranges
    RooRealVar var_IM("IM","IM", h_data->GetXaxis()->GetXmin(), h_data->GetXaxis()->GetXmax(), "MeV");
    var_IM.setBins(10000);
    var_IM.setRange("full",var_IM.getMin(),var_IM.getMax());

    // load data to be fitted
    RooDataHist h_roo_data("h_roo_data","dataset",var_IM,h_data);

    // build shifted mc lineshape
    RooRealVar var_IM_shift("var_IM_shift", "shift in IM", 3.0, -10.0, 10.0);
    RooProduct var_IM_shift_invert("var_IM_shift_invert","shifted IM",RooArgSet(var_IM_shift, RooConst(-1.0)));
    RooAddition var_IM_shifted("var_IM_shifted","shifted IM",RooArgSet(var_IM,var_IM_shift_invert));
    RooDataHist h_roo_mc("h_roo_mc","MC lineshape", var_IM, h_mc);
    RooHistPdf pdf_mc_lineshape("pdf_mc_lineshape","MC lineshape as PDF", var_IM_shifted, var_IM, h_roo_mc, 4);

    // build detector resolution smearing
    RooRealVar  var_gauss_sigma_1("gauss_sigma_1","width of gaussian 1",  4.0, 0.0, 100.0);
    RooRealVar  var_gauss_sigma_2("gauss_sigma_2","width of gaussian 2", 10.0, 0.0, 100.0);
    RooRealVar  var_gauss_fraction("gauss_fraction","fraction for gaussian 1/2", 0.9, 0.0, 1.0);

//    const double width = var_IM.getMax()-var_IM.getMin();
//    RooRealVar var_kernel("var_kernel","var_kernel", -width/2, +width/2 , "MeV");
//    var_kernel.setBins(10000);
    RooGaussian pdf_gaussian1("pdf_gaussian_1","Gaussian 1", var_IM, RooConst(0.0), var_gauss_sigma_1);
    RooGaussian pdf_gaussian2("pdf_gaussian_2","Gaussian 2", var_IM, RooConst(0.0), var_gauss_sigma_2);

    // double gaussian
//    RooAddPdf pdf_smearing("pdf_smearing","Double Gaussian",RooArgList(pdf_gaussian1,pdf_gaussian2),var_gauss_fraction);

    // single gaussian
    RooGaussian pdf_smearing("pdf_smearing","Single Gaussian", var_IM, RooConst(0.0), var_gauss_sigma_1);

    // Crystal Ball shap
//    RooRealVar var_
//    RooCBShape pdf_smearing("pdf_smearing", "CBShape", )

    // build signal as convolution, note that the gaussian must be the second PDF (see documentation)
    RooFFTConvPdf pdf_signal("pdf_signal","MC_lineshape (X) gauss", var_IM, pdf_mc_lineshape, pdf_smearing);


    // build signal as simple gaussian (not working actually)

//    RooRealVar  var_gauss_mean("gauss_mean","mean of gaussian", 958.0, 900.0, 1000.0);
//    RooGaussian pdf_signal("pdf_signal","Gaussian signal", var_IM, var_gauss_mean, var_gauss_sigma);


    // build background (chebychev or argus?)

//    const int polOrder = 6;
//    std::vector<std::unique_ptr<RooRealVar>> bkg_params; // RooRealVar cannot be copied, so create them on heap
//    RooArgSet roo_bkg_params;
//    for(int p=0;p<polOrder;p++) {
//        bkg_params.emplace_back(std_ext::make_unique<RooRealVar>((
//                                    "p_"+to_string(p)).c_str(), ("Bkg Par "+to_string(p)).c_str(), 0.0, -10.0, 10.0)
//                                );
//        roo_bkg_params.add(*bkg_params.back());
//    }
//    RooChebychev pdf_background("chebychev","Polynomial background",var_IM,roo_bkg_params);
//    var_IM.setRange("bkg_l", var_IM.getMin(), 930);
//    var_IM.setRange("bkg_r", 990, 1000);
//    pdf_background.fitTo(h_roo_data, Range("bkg_l,bkg_r"), Extended()); // using Range(..., ...) does not work here (bug in RooFit, sigh)

    RooRealVar argus_cutoff("argus_cutoff","argus pos param", 1004, 1000, 1050);
    RooRealVar argus_shape("argus_shape","argus shape param", -15, -50.0, 0.0);
    RooRealVar argus_p("argus_p","argus p param", 3.4, 0.0, 4.0);
    RooArgusBG pdf_background("argus","bkg argus",var_IM,argus_cutoff,argus_shape,argus_p);

    // build sum
    RooRealVar nsig("nsig","#signal events", 5e4, 0, 1e5);
    RooRealVar nbkg("nbkg","#background events", 8e4, 0, 1e5);
    RooAddPdf pdf_sum("pdf_sum","total sum",RooArgList(pdf_signal,pdf_background),RooArgList(nsig,nbkg));

    // do some pre-fitting to obtain better starting values, make sure function is non-zero in range
//    var_IM.setRange("nonzero",var_IM.getMin(), 1000.0);
//    pdf_sum.chi2FitTo(h_roo_data, Range("nonzero"), PrintLevel(-1), Optimize(false)); // using Range(..., ...) does not work here (bug in RooFit, sigh)


//    pdf_signal.Print("v");

//    RooTrace::dump() ;

//    RooPlot* frame_test = var_kernel.frame();
//    pdf_sum.plotOn(frame_test);
//    frame_test->Draw();

//    pdf_sum.Print("v");


    // do the actual maximum likelihood fit
    auto fr_data = pdf_sum.fitTo(h_roo_data, Extended(), SumW2Error(kTRUE), Range("full"), Save());
    const auto numParams = fr_data->floatParsFinal().getSize();

    // draw output, won't be shown in batch mode
    RooPlot* frame_data = var_IM.frame();
    h_roo_data.plotOn(frame_data);
//    pdf_sum.plotOn(frame, LineColor(kRed), VisualizeError(*fr));
    pdf_sum.plotOn(frame_data, LineColor(kRed));
    const auto chi2ndf = frame_data->chiSquare(numParams);
    auto pdf_sum_tf = pdf_sum.asTF(var_IM);
    const auto peak_pos = pdf_sum_tf->GetMaximumX(signal_region.Start(), signal_region.Stop());
    auto hresid = frame_data->residHist();

    pdf_sum.plotOn(frame_data, Components(pdf_background), LineColor(kBlue));
    pdf_sum.plotOn(frame_data, Components(pdf_signal), LineColor(kGreen));

    frame_data->Draw();

    const auto text_x = peak_pos+signal_region.Length()*0.1;
    const auto text_y = h_data->GetMaximum();
    auto text = new TPaveText(text_x, text_y-h_data->GetMaximum()*0.3, text_x+signal_region.Length()*1.0, text_y);
    text->SetBorderSize(0);
    text->SetFillColor(kWhite);
    text->AddText(static_cast<string>(std_ext::formatter() << "N_{sig} = " << nsig.getVal() << " #pm " << nsig.getError()).c_str());
    text->AddText(static_cast<string>(std_ext::formatter() << "#chi^{2}_{red} = " << chi2ndf).c_str());
    text->Draw();

//    RooPlot* frame2 = var_IM.frame(Title("Residual Distribution")) ;
//    frame2->addPlotable(hresid,"P");
//    new TCanvas();
//    frame2->Draw();

    fr_data->Print("v");

    LOG(INFO) << "peakPos=" << peak_pos;
    LOG(INFO) << "numParams=" << numParams << " chi2ndf=" << chi2ndf;
    LOG(INFO) << "residuals_integral/perbin="
              << hresid->Integral(hresid->GetXaxis()->FindBin(signal_region.Start()),
                                  hresid->GetXaxis()->FindBin(signal_region.Stop()))/hresid->getNominalBinWidth();


    // do fit on MC input



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
