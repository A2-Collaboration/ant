#include "base/Logger.h"

#include "tclap/CmdLine.h"
#include "tclap/ValuesConstraintExtra.h"

#include "base/interval.h"
#include "base/WrapTFile.h"
#include "base/std_ext/string.h"
#include "base/std_ext/system.h"
#include "base/std_ext/memory.h"
#include "base/std_ext/math.h"
#include "base/ParticleType.h"
#include "base/TH_ext.h"

#include "expconfig/ExpConfig.h"
#include "base/Detector_t.h"


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
//#include "RooRealVar.h"
//#include "RooGaussian.h"
//#include "RooArgusBG.h"
//#include "RooAddPdf.h"
//#include "RooDataSet.h"
//#include "RooHistPdf.h"
//#include "RooHist.h"
//#include "RooPlot.h"
//#include "RooDataHist.h"
//#include "RooAddition.h"
//#include "RooProduct.h"
//#include "RooChebychev.h"
//#include "RooConstVar.h"
//#include "RooDerivative.h"
//#include "RooFFTConvPdf.h"
//#include "RooChi2Var.h"
//#include "RooMinuit.h"
//#include "RooFitResult.h"


using namespace ant;
using namespace std;
//using namespace RooFit;

int main(int argc, char** argv) {
    SetupLogger();

    TCLAP::CmdLine cmd("SinglePi0_fit", ' ', "0.1");

    auto cmd_verbose       = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"int");

    auto cmd_data          = cmd.add<TCLAP::ValueArg<string>>("","data","Data input from  singlePi0-Plot",true,"","rootfile");

    auto cmd_lumi          = cmd.add<TCLAP::ValueArg<string>>("","lumi","path to luminosity-class output if seperate file from data input",true,"","rootfile");

    auto cmd_eff           = cmd.add<TCLAP::ValueArg<string>>("","eff","MC signal input",true,"","rootfile");

    auto cmd_histpath      = cmd.add<TCLAP::ValueArg<string>>("","histpath","Path for hists (determines cutstr)",false,
                                                              "dicardedEk<20/EMB_prob>0.05/Pi0PIDVeto==0","path");

    auto cmd_histname      = cmd.add<TCLAP::ValueArg<string>>("","histname","Name of hist",false,"recon_fit","name");

    auto cmd_histluminame  = cmd.add<TCLAP::ValueArg<string>>("","histlumi","Name of hist",false,"intlumicor","name");

    auto cmd_histreconame  = cmd.add<TCLAP::ValueArg<string>>("","histreco","Name of hist",false,"recon_fit","name");
    auto cmd_histseenname  = cmd.add<TCLAP::ValueArg<string>>("","histseen","Name of hist",false,"seenMCcosTheta","name");

    TCLAP::ValuesConstraintExtra<decltype(ExpConfig::Setup::GetNames())> allowedsetupnames(ExpConfig::Setup::GetNames());
    auto cmd_setup  = cmd.add<TCLAP::ValueArg<string>>("s","setup","Choose setup by name",true,"", &allowedsetupnames);


    auto cmd_batchmode = cmd.add<TCLAP::MultiSwitchArg>("b", "batch","Run in batch mode (no ROOT shell afterwards)", false);
    auto cmd_mchmode   = cmd.add<TCLAP::MultiSwitchArg>("",  "mc",   "run on mc",                                    false);

    auto cmd_output    = cmd.add<TCLAP::ValueArg<string>>("o","output","Output file",false,"","filename");


    cmd.parse(argc, argv);
    if(cmd_verbose->isSet()) {
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());
    }
    auto loadHist = [](const WrapTFileInput& input, const string& histpath)
    {
        TH1* hist = nullptr;
        if (!input.GetObject(histpath,hist))
            throw runtime_error(std_ext::formatter() << "Cannot find " << histpath);
        return hist;
    };

    ExpConfig::Setup::SetByName(cmd_setup->getValue());
    auto Tagger = ExpConfig::Setup::GetDetector<TaggerDetector_t>();

    const string cosThetaLabel = "cos(#theta_{cm})";
    const string taggerLabel   = "tagger channel";
    const string xsecLabel     = "#frac{d#sigma}{d#Omega} [#mub/sr]";

    const pair<double,double> userRangeTheta({-0.9,0.9});

    const string plotterPath = "singlePi0_Plot/";
    WrapTFileInput input_data(cmd_data->getValue());
    WrapTFileInput input_eff(cmd_eff->getValue());

    const string dataSource = cmd_mchmode->isSet() ? "/h/Sum_MC/" : "/h/data/";
    auto h_data = dynamic_cast<TH2D*>(loadHist(input_data,
                                               plotterPath + cmd_histpath->getValue() + dataSource +cmd_histname->getValue()));
    auto h_rec  = dynamic_cast<TH2D*>(loadHist(input_eff,
                                               plotterPath + cmd_histpath->getValue() + "/h/Sig/" + cmd_histreconame->getValue()));
    auto h_seen = dynamic_cast<TH2D*>(loadHist(input_eff,
                                               plotterPath + cmd_histseenname->getValue()));


    const string fluxPath    = "PhotonFlux/";
    TH1D* h_lumi = nullptr;
    if (cmd_mchmode->isSet())
    {
        LOG(INFO) << "Explicitly set to 'mc', using L == 1" ;
    } else
    {
        WrapTFileInput input_lumi(cmd_lumi->getValue());
        h_lumi = dynamic_cast<TH1D*>(loadHist(input_lumi,
                                              fluxPath + cmd_histluminame->getValue()));
    }
    const auto nChannels       = h_data->GetNbinsX();
    if (nChannels != static_cast<int>(Tagger->GetNChannels()))
    {
        LOG(ERROR) << "hitograms don't match with provided setup";
        exit(1);
    }
    const BinSettings taggBins = BinSettings(nChannels);
    const auto cosThetaBins    = TH_ext::getBins(h_data->GetYaxis());
    const auto DeltaOmega      = cosThetaBins.Length() * 2 * M_PI / cosThetaBins.Bins();


    argc=0; // prevent TRint to parse any cmdline
    TRint app("SinglePi0_fit",&argc,argv,nullptr,0,true);

    unique_ptr<WrapTFileOutput> masterFile;
    if(cmd_output->isSet()) {
        // cd into masterFile upon creation
        masterFile = std_ext::make_unique<WrapTFileOutput>(cmd_output->getValue(), true);
    }

    analysis::HistogramFactory histfac("singlePi0_fits");
    auto histEff = histfac.makeTH2D("eff",
                                    taggerLabel,cosThetaLabel,
                                    taggBins,cosThetaBins,"eff",true);
    histEff->Add(h_rec);
    histEff->Divide(h_seen);

    auto lumi2d  = histfac.makeTH2D("luminosity",
                                    taggerLabel,cosThetaLabel,
                                    taggBins,cosThetaBins,
                                    "lumi2d",true);
    for (int tagBin = 1 ; tagBin <= nChannels ; ++tagBin)
    {
        for (int thetaBin = 1 ; thetaBin <= static_cast<int>(cosThetaBins.Bins()) ; ++thetaBin)
        {
            if (h_lumi)
            {
                lumi2d->SetBinContent(tagBin,thetaBin,h_lumi->GetBinContent(tagBin));
                lumi2d->SetBinError(tagBin,thetaBin,h_lumi->GetBinError(tagBin));
            } else
            {
                lumi2d->SetBinContent(tagBin,thetaBin,1.0);
                lumi2d->SetBinError(tagBin,thetaBin,0.0);
            }
        }
    }

    auto sigma2d = histfac.makeTH2D("cross sections",
                                    taggerLabel,cosThetaLabel,
                                    taggBins,cosThetaBins,
                                    "sigma2d",true);
    sigma2d->Add(h_data);
    sigma2d->Divide(histEff);
    sigma2d->Divide(lumi2d);
    sigma2d->Scale(1./DeltaOmega);




    auto applyCosmetics = [&userRangeTheta,&xsecLabel] (TH1D* hist, const string& title, const bool isCosTheta = true)
    {
        hist->SetTitle(title.c_str());
        if (isCosTheta)
            hist->GetXaxis()->SetRangeUser(userRangeTheta.first,userRangeTheta.second);
        hist->GetYaxis()->SetTitle(xsecLabel.c_str());
    };

    vector<TH1D*> histChannels(nChannels);
    for (auto ch = 0 ; ch < nChannels ; ++ch)
    {
        const string hname = std_ext::formatter() << "ch" << ch;

        const auto energy_interval = IntervalD::CenterWidth(Tagger->GetPhotonEnergy(ch),Tagger->GetPhotonEnergyWidth(ch));


        histChannels[ch] = sigma2d->ProjectionY(hname.c_str(),ch+1,ch+1);
        applyCosmetics(histChannels[ch],std_ext::formatter() << "Differential cross section: E_{#gamma} in " << energy_interval << " MeV");
    }

    auto histsigma_Theta = sigma2d->ProjectionY("SigmaTheta");
    histsigma_Theta->Scale(1. / nChannels);
    applyCosmetics(histsigma_Theta,
                   std_ext::formatter() << "Differential cross section for E_{#gamma} in "
                                        << IntervalD(Tagger->GetPhotonEnergy(nChannels-1) - Tagger->GetPhotonEnergyWidth(nChannels-1) / 2,
                                                     Tagger->GetPhotonEnergy(0) + Tagger->GetPhotonEnergyWidth(0) / 2)
                                        << " MeV");
    auto histsigma_E = sigma2d->ProjectionX("sigmaE",4,28,"widthq");
    applyCosmetics(histsigma_E,
                   "Total cross sections",
                   false);


    for (auto h: histChannels)
        h->SetBit(TH1::kIsAverage);
    histsigma_Theta->SetBit(TH1::kIsAverage);
    histsigma_E->SetBit(TH1::kIsAverage);


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
