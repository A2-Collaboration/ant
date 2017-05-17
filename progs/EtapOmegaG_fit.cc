#include "base/Logger.h"

#include "tclap/CmdLine.h"
#include "base/interval.h"
#include "base/WrapTFile.h"
#include "base/std_ext/string.h"
#include "base/std_ext/system.h"
#include "base/std_ext/memory.h"
#include "base/std_ext/math.h"
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

using namespace ant;
using namespace std;
using namespace RooFit;

struct fit_params_t {
    interval<double> signal_region{920, 990};
    int nSamplingBins{10000};
    int interpOrder{4};

    TH1D* h_mc = nullptr;
    TH1D* h_data = nullptr;

};

struct fit_return_t : ant::root_drawable_traits {

    fit_params_t p;

    RooFitResult* fitresult = nullptr;
    double chi2ndf = std_ext::NaN;
    double peakpos = std_ext::NaN;

    int numParams() {
        return fitresult->floatParsFinal().getSize();
    }

    double residualSignalIntegral() {
        auto& h = residual;
        return h->Integral(h->GetXaxis()->FindBin(p.signal_region.Start()),
                           h->GetXaxis()->FindBin(p.signal_region.Stop()))
                /h->getNominalBinWidth();
    }

    RooPlot* fitplot = nullptr;
    TPaveText* lbl = nullptr;
    RooHist* residual = nullptr;

    virtual void Draw(const string& option) const override
    {
        auto& pars = fitresult->floatParsFinal();
        auto nsig = (RooRealVar*)pars.at(pars.index("nsig"));

        (void)option;

        const auto text_x = peakpos+p.signal_region.Length()*0.1;
        const auto text_y = p.h_data->GetMaximum();
        auto lbl = new TPaveText(text_x, text_y-p.h_data->GetMaximum()*0.3, text_x+p.signal_region.Length()*1.0, text_y);
        lbl->SetBorderSize(0);
        lbl->SetFillColor(kWhite);
        lbl->AddText(static_cast<string>(std_ext::formatter() << "N_{sig} = " << nsig->getVal() << " #pm " << nsig->getError()).c_str());
        lbl->AddText(static_cast<string>(std_ext::formatter() << "#chi^{2}_{red} = " << chi2ndf).c_str());

        fitplot->Draw();
        lbl->Draw();
    }
};

fit_return_t doFit(const fit_params_t& p) {

    fit_return_t r;
    r.p = p; // remember params

    // define observable and ranges
    RooRealVar var_IM("IM","IM", p.h_data->GetXaxis()->GetXmin(), p.h_data->GetXaxis()->GetXmax(), "MeV");
    var_IM.setBins(p.nSamplingBins);
    var_IM.setRange("full",var_IM.getMin(),var_IM.getMax());

    // load data to be fitted
    RooDataHist h_roo_data("h_roo_data","dataset",var_IM,p.h_data);

    // build shifted mc lineshape
    RooRealVar var_IM_shift("var_IM_shift", "shift in IM", 3.0, -10.0, 10.0);
    RooProduct var_IM_shift_invert("var_IM_shift_invert","shifted IM",RooArgSet(var_IM_shift, RooConst(-1.0)));
    RooAddition var_IM_shifted("var_IM_shifted","shifted IM",RooArgSet(var_IM,var_IM_shift_invert));
    RooDataHist h_roo_mc("h_roo_mc","MC lineshape", var_IM, p.h_mc);
    RooHistPdf pdf_mc_lineshape("pdf_mc_lineshape","MC lineshape as PDF", var_IM_shifted, var_IM, h_roo_mc, p.interpOrder);

    // build detector resolution smearing

    RooRealVar  var_sigma("sigma","detector resolution",  4.0, 0.0, 100.0);
    RooGaussian pdf_smearing("pdf_smearing","Single Gaussian", var_IM, RooConst(0.0), var_sigma);

    // build signal as convolution, note that the gaussian must be the second PDF (see documentation)
    RooFFTConvPdf pdf_signal("pdf_signal","MC_lineshape (X) gauss", var_IM, pdf_mc_lineshape, pdf_smearing);

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

    RooRealVar argus_cutoff("argus_cutoff","argus pos param", 1033, 1000, 1050);
    RooRealVar argus_shape("argus_shape","argus shape param", -15, -50.0, 0.0);
    RooRealVar argus_p("argus_p","argus p param", 3.4, 0.0, 4.0);
    RooArgusBG pdf_background("argus","bkg argus",var_IM,argus_cutoff,argus_shape,argus_p);

    // build sum
    RooRealVar nsig("nsig","#signal events", 7e4, 0, 1e6);
    RooRealVar nbkg("nbkg","#background events", 1e5, 0, 1e6);
    RooAddPdf pdf_sum("pdf_sum","total sum",RooArgList(pdf_signal,pdf_background),RooArgList(nsig,nbkg));

    // do some pre-fitting to obtain better starting values, make sure function is non-zero in range
    //    var_IM.setRange("nonzero",var_IM.getMin(), 1000.0);
    //    pdf_sum.chi2FitTo(h_roo_data, Range("nonzero"), PrintLevel(-1), Optimize(false)); // using Range(..., ...) does not work here (bug in RooFit, sigh)


    // do the actual maximum likelihood fit
    // use , Optimize(false), Strategy(2) for double gaussian...?!
    r.fitresult = pdf_sum.fitTo(h_roo_data, Extended(), SumW2Error(kTRUE), Range("full"), Save(), PrintLevel(-1));

    // draw output and remember pointer
    r.fitplot = var_IM.frame();

    h_roo_data.plotOn(r.fitplot);
    //    pdf_sum.plotOn(frame, LineColor(kRed), VisualizeError(*fr));

    // need to figure out chi2nds and stuff after plotting data and finally fitted pdf_sum
    // also the residHist must be created here (and rememebered for later use)
    pdf_sum.plotOn(r.fitplot, LineColor(kRed));
    r.chi2ndf = r.fitplot->chiSquare(r.numParams());
    auto pdf_sum_tf = pdf_sum.asTF(var_IM);
    r.peakpos = pdf_sum_tf->GetMaximumX(p.signal_region.Start(), p.signal_region.Stop());
    r.residual = r.fitplot->residHist();

    pdf_sum.plotOn(r.fitplot, Components(pdf_background), LineColor(kBlue));
    pdf_sum.plotOn(r.fitplot, Components(pdf_signal), LineColor(kGreen));

    return r;
}

int main(int argc, char** argv) {
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

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

    fit_params_t p;

    WrapTFileInput input_data(cmd_data->getValue()); // keep it open
    {
        const string histpath = cmd_histpath->getValue()+"/h/Data/"+cmd_histname->getValue();
        if(!input_data.GetObject(histpath, p.h_data)) {
            LOG(ERROR) << "Cannot find " << histpath;
            return EXIT_FAILURE;
        }
    }

    WrapTFileInput input_mc(cmd_mc->getValue()); // keep it open
    {
        const string histpath = cmd_histpath->getValue()+"/h/Sum_MC/"+cmd_histname->getValue();
        if(!input_mc.GetObject(histpath, p.h_mc)) {
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

    auto r = doFit(p);


    //    RooPlot* frame2 = var_IM.frame(Title("Residual Distribution")) ;
    //    frame2->addPlotable(hresid,"P");
    //    new TCanvas();
    //    frame2->Draw();

//    r.fitresult->Print("v");

//    LOG(INFO) << "peakPos=" << r.peakpos;
//    LOG(INFO) << "numParams=" << r.numParams() << " chi2ndf=" << r.chi2ndf;
//    LOG(INFO) << "residuals_integral/perbin=" << r.residualSignalIntegral();

    ant::canvas("EtapOmegaG_fit") << r << endc;

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
