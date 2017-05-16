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


    RooRealVar var_IM("IM","IM", h_data->GetXaxis()->GetXmin(), h_data->GetXaxis()->GetXmax(), "MeV");
    RooDataHist data("h_data","dataset",var_IM,h_data);

//    // --- Build Gaussian signal PDF ---
//    RooRealVar sigmean("sigmean","B^{#pm} mass",5.28,5.20,5.30);
//    RooRealVar sigwidth("sigwidth","B^{#pm} width",0.0027,0.001,1.);
//    RooGaussian gauss("gauss","gaussian PDF",mes,sigmean,sigwidth);

//    // --- Build Argus background PDF ---
//    RooRealVar argpar1("argpar1","argus parameter 1",5.291,0.0,10.0);
//    RooRealVar argpar2("argpar2","argus shape parameter",-20.0,-100.,-1.);
//    RooArgusBG argus("argus","Argus PDF",mes,argpar1,argpar2);

//    // --- Construct signal+background PDF ---
//    RooRealVar nsig("nsig","#signal events",200,0.,10000);
//    RooRealVar nbkg("nbkg","#background events",800,0.,10000);
//    RooAddPdf sum("sum","g+a",RooArgList(gauss,argus),RooArgList(nsig,nbkg));

//    // --- Generate a toyMC sample from composite PDF ---
//    RooDataSet *data = sum.generate(mes,2000);
//    // --- Perform extended ML fit of composite PDF to toy data ---
//    sum.fitTo(*data,Extended());


    if(!cmd_batchmode->isSet()) {
        if(!std_ext::system::isInteractive()) {
            LOG(INFO) << "No TTY attached. Not starting ROOT shell.";
        }
        else {

            argc=0; // prevent TRint to parse any cmdline
            TRint app("EtapOmegaG_plot",&argc,argv,nullptr,0,true);

            if(masterFile)
                LOG(INFO) << "Close ROOT properly to write data to disk.";

            // --- Plot toy data and composite PDF overlaid ---
            RooPlot* frame = var_IM.frame();
            data.plotOn(frame);
            frame->Draw();

            app.Run(kTRUE); // really important to return...
            if(masterFile)
                LOG(INFO) << "Writing output file...";
            masterFile = nullptr;   // and to destroy the master WrapTFile before TRint is destroyed
        }
    }

    return 0;
}
