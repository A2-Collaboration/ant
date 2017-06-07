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
    auto cmd_histname      = cmd.add<TCLAP::ValueArg<string>>("","histname","Name of hist",false,"recon_cor","name");

    auto cmd_histluminame  = cmd.add<TCLAP::ValueArg<string>>("","histlumi","Name of hist",false,"intlumi","name");

    auto cmd_histreconame  = cmd.add<TCLAP::ValueArg<string>>("","histreco","Name of hist",false,"effrecon_pi0","name");
    auto cmd_histseenname  = cmd.add<TCLAP::ValueArg<string>>("","histseen","Name of hist",false,"seenMCcosTheta","name");



    auto cmd_batchmode = cmd.add<TCLAP::MultiSwitchArg>("b","batch","Run in batch mode (no ROOT shell afterwards)",false);
    auto cmd_output    = cmd.add<TCLAP::ValueArg<string>>("o","output","Output file",false,"","filename");


    cmd.parse(argc, argv);
    if(cmd_verbose->isSet()) {
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());
    }


    TH2D* h_data = nullptr;
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


    TH2D* h_rec= nullptr;
    TH2D* h_seen = nullptr;
    WrapTFileInput input_eff(cmd_eff->getValue()); // keep it open
    {
        const string histpath_rec  = cmd_histpath->getValue() + "/h/Sig/" + cmd_histreconame->getValue();
        const string histpath_seen = "singlePi0_Plot/" + cmd_histseenname->getValue();

        if(!input_eff.GetObject(histpath_rec, h_rec)) {
            LOG(ERROR) << "Cannot find " << histpath_rec;
            return EXIT_FAILURE;
        }
        if(!input_eff.GetObject(histpath_seen, h_seen)) {
            LOG(ERROR) << "Cannot find " << histpath_seen;
            return EXIT_FAILURE;
        }
    }

    const auto nChannels       = h_data->GetNbinsX();
    const auto cosThetaBinning = TH_ext::getBins(h_data->GetYaxis());
    const auto DeltaOmega      = cosThetaBinning.Length() * 2 * M_PI / cosThetaBinning.Bins();



    argc=0; // prevent TRint to parse any cmdline
    TRint app("OmegaEtaG_fit",&argc,argv,nullptr,0,true);

    unique_ptr<WrapTFileOutput> masterFile;
    if(cmd_output->isSet()) {
        // cd into masterFile upon creation
        masterFile = std_ext::make_unique<WrapTFileOutput>(cmd_output->getValue(), true);
    }
    analysis::HistogramFactory histfac("singlePi0_fits");

    auto histEff = histfac.makeTH2D("eff",
                                    "tagger channel","cos(#theta)",
                                    BinSettings(nChannels),cosThetaBinning,"eff",true);
    histEff->Add(h_rec);
    histEff->Divide(h_seen);

    auto histAll = histfac.makeTH1D("all channels","cos(#theta_{coms})","#frac{d#sigma}{d#Omega} [#mub/sr]",
                                    cosThetaBinning,"all",true);

    vector<TH1D*> histChannels(nChannels);
    for (auto ch = 0 ; ch < nChannels ; ++ch)
    {
        const string hname = std_ext::formatter() << "ch" << ch;

        histChannels[ch] = histfac.makeTH1D(hname,"cos(#theta_{coms})","#frac{d#sigma}{d#Omega} [#mub/sr]",
                                            cosThetaBinning,hname,true);
        histChannels[ch]->Add(h_data->ProjectionY("py_d",ch+1,ch+1));
        histAll->Add(histChannels[ch]);
        histChannels[ch]->Divide(histEff->ProjectionY("py_e",ch+1,ch+1));
        histChannels.at(ch)->Scale(1./(DeltaOmega * h_lumi->GetBinContent(ch+1)));
    }
    histAll->Divide(histEff->ProjectionY());
    histAll->Scale(1./(DeltaOmega * h_lumi->Integral()));

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
