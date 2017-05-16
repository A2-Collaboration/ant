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
#include "RooChebychev.h"

using namespace ant;
using namespace std;
using namespace RooFit;

int main(int argc, char** argv) {
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

    unique_ptr<WrapTFileOutput> masterFile;
    if(cmd_output->isSet()) {
        // cd into masterFile upon creation
        masterFile = std_ext::make_unique<WrapTFileOutput>(cmd_output->getValue(), true);
    }




    // define observable and ranges
    RooRealVar var_IM("IM","IM", h_data->GetXaxis()->GetXmin(), h_data->GetXaxis()->GetXmax(), "MeV");
    var_IM.setRange("bkg_l", var_IM.getMin(), 930);
    var_IM.setRange("bkg_r", 990, var_IM.getMax());
    var_IM.setRange("nonzero", var_IM.getMin(), 1000);


    // load data to be fitted
    RooDataHist h_roo_data("h_roo_data","dataset",var_IM,h_data);

    // build shifted mc lineshape
    RooRealVar var_IM_shift("var_IM_shift", "shift in IM", -3.0, -10, 10);
    RooAddition var_IM_shifted("var_IM_shifted","shifted IM",RooArgSet(var_IM,var_IM_shift));
    RooDataHist h_roo_mc("h_roo_mc","MC lineshape", var_IM, h_mc);
    RooHistPdf pdf_mc_lineshape("pdf_mc_lineshape","MC lineshape as PDF", var_IM_shifted, var_IM, h_roo_mc, 4);

    // build background
//    const int polOrder = 5;
//    std::vector<std::unique_ptr<RooRealVar>> bkg_params; // RooRealVar cannot be copied, so create them on heap
//    RooArgSet roo_bkg_params;
//    for(int p=0;p<polOrder;p++) {
//        bkg_params.emplace_back(std_ext::make_unique<RooRealVar>((
//                                    "p_"+to_string(p)).c_str(), ("Bkg Par "+to_string(p)).c_str(), 0, 0, 1)
//                                );
//        roo_bkg_params.add(*bkg_params.back());
//    }
//    RooChebychev pdf_background("pdf_background","Polynomial background",var_IM,roo_bkg_params);

    RooRealVar argus_pos("argus_pos","argus pos param", 1010, 1000, 1050);
    RooRealVar argus_shape("argus_shape","argus shape param", -5.0, -10.0, 0.0);
    RooRealVar argus_p("argus_p","argus p param", 1.2, 0.0, 2.0);
    RooArgusBG pdf_background("argus","bkg argus",var_IM,argus_pos,argus_shape,argus_p);


    // build sum
    RooRealVar nsig("nsig","#signal events", 600, 0, 100000);
    RooRealVar nbkg("nbkg","#background events", 1000, 0, 100000);
    RooAddPdf pdf_sum("pdf_sum","total sum",RooArgList(pdf_mc_lineshape,pdf_background),RooArgList(nsig,nbkg));

    if(!cmd_batchmode->isSet()) {
        if(!std_ext::system::isInteractive()) {
            LOG(INFO) << "No TTY attached. Not starting ROOT shell.";
        }
        else {

            argc=0; // prevent TRint to parse any cmdline
            TRint app("EtapOmegaG_plot",&argc,argv,nullptr,0,true);

            if(masterFile)
                LOG(INFO) << "Close ROOT properly to write data to disk.";

            // do some fitting
            pdf_sum.chi2FitTo(h_roo_data, Range("nonzero")); // using Range(..., ...) does not work here!

            RooPlot* frame = var_IM.frame();
            h_roo_data.plotOn(frame);
            pdf_sum.plotOn(frame);
            pdf_sum.plotOn(frame, Components(pdf_background), LineStyle(kDashed)) ;
            frame->Draw();

            app.Run(kTRUE); // really important to return...
            if(masterFile)
                LOG(INFO) << "Writing output file...";
            masterFile = nullptr;   // and to destroy the master WrapTFile before TRint is destroyed
        }
    }

    return 0;
}
