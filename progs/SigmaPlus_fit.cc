#include "base/Logger.h"

#include "tclap/CmdLine.h"
#include "tclap/ValuesConstraintExtra.h"
#include "base/WrapTFile.h"
#include "base/ParticleType.h"
#include "analysis/utils/ValError.h"
#include "expconfig/ExpConfig.h"
#include "analysis/plot/HistogramFactory.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TGraphErrors.h"

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
#include "RooCBShape.h"
#include "RooProdPdf.h"
#include "RooExtendPdf.h"



using namespace ant;
using namespace std;
using namespace RooFit;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::utils;

unsigned globalBin = 0;


ValError MakeDataPoint(const RooRealVar& value)
{
    return ValError(value.getValV(),value.getError());
}

struct CrossSectionDataPoint
{
    double Egamma;
    ValError Data;
    string Unit;

    CrossSectionDataPoint(const double eGamma, const ValError data, const string& unit = "#mub"):
        Egamma(eGamma), Data(data), Unit(unit) {}
    CrossSectionDataPoint(const double eGamma, const RooRealVar& value, const string& unit = "#mub"):
        Egamma(eGamma), Data(MakeDataPoint(value)), Unit(unit) {}

    double W() const {return sqrt( sqr(Egamma + ParticleTypeDatabase::Proton.Mass()) - sqr(Egamma));}
};

ValError fitHist (TH2D* histDalitz) {

  RooRealVar mppi0("mppi0","m(p #pi^{0})",1150, 1230) ;
  RooRealVar mpi0pi0("mpi0pi0","m(#pi^{0} #pi^{0})", 420, 560) ;

  RooDataHist dh ("dh", "dh", RooArgSet(mppi0, mpi0pi0), histDalitz);


  RooRealVar mean_Sigma("mean_Sigma","m_{#Sigma}",1188, 1170, 1220);
  RooRealVar sigma_Sigma("sigma_Sigma","#sigma_{#Sigma}",10, 0, 50);
  RooRealVar alpha_Sigma("alpha_Sigma","#alpha_{#Sigma}",2);
  RooRealVar n_Sigma("n_Sigma","n_Sigma",1);
  RooCBShape *theSigma = new RooCBShape("theSigma","crystal ball PDF",mppi0,mean_Sigma,sigma_Sigma,alpha_Sigma,n_Sigma);

  RooRealVar p0_Sigma ("p0_Sigma", "p0_Sigma", 0, -1, 1);
  RooRealVar p1_Sigma ("p1_Sigma", "p1_Sigma", 0, -1, 1);
  RooRealVar p2_Sigma ("p2_Sigma", "p2_Sigma", 0, -1, 1);
  RooRealVar p3_Sigma ("p3_Sigma", "p3_Sigma", 0);

  RooRealVar frac_peak ("frac_peak", "frac_peak", 0.5, 0, 1);

  RooChebychev *bkg_Sigma = new RooChebychev ("bkg_Sigma", "bkg_Sigma",
                          mppi0, RooArgList(p0_Sigma, p1_Sigma, p2_Sigma, p3_Sigma));

  //  RooAddPdf *SigmaPdf = new RooAddPdf ("SigmaPdf", "SigmaPdf",
  //				       RooArgList (*theSigma, *bkg_Sigma),
  //				       RooArgList (frac_peak));

  RooRealVar mean_K0S("mean_K0S","m_{K_{S}}",490, 460, 500);
  RooRealVar sigma_K0S("sigma_K0S","#sigma_{K_{S}}",12, 0, 50);
  RooRealVar alpha_K0S("alpha_K0S","#alpha_{K_{S}}",0.97);
  RooRealVar n_K0S("n_K0S","n_K0S",10);
  RooCBShape *theK0S = new RooCBShape("theK0S","crystal ball PDF for K0S",mpi0pi0,mean_K0S,sigma_K0S,alpha_K0S,n_K0S);

  RooRealVar p0_K0S ("p0_K0S", "p0_K0S", 0, -1, 1);
  RooRealVar p1_K0S ("p1_K0S", "p1_K0S", 0, -1, 1);
  RooRealVar p2_K0S ("p2_K0S", "p2_K0S", 0, -1, 1);

  RooRealVar frac_peak_K0S ("frac_peak_K0S", "frac_peak_K0S", 0.5, 0, 1);

  RooChebychev *bkg_K0S = new RooChebychev ("bkg_K0S", "bkg_K0S",
                        mpi0pi0, RooArgList(p0_K0S, p1_K0S, p2_K0S));


  //RooAddPdf *K0SPdf = new RooAddPdf ("K0SPdf", "K0SPdf",
  //				     RooArgList (*theK0S, *bkg_K0S),
  //				     RooArgList (frac_peak_K0S));


  //  SigmaPdf->fitTo (dh);
  //  K0SPdf  ->fitTo (dh);

  RooProdPdf *sigPdf = new RooProdPdf ("sigPdf", "signal 2D PDF", RooArgList (*theSigma, *theK0S));
  RooProdPdf *bkgPdf = new RooProdPdf ("bkgPdf", "background 2D PDF", RooArgList (*bkg_Sigma, *bkg_K0S));

  RooRealVar nSig ("nSig", "nSig", 10000,  0, 10000000);
  RooRealVar nBkg ("nBkg", "nBkg", 100000, 0, 10000000);

  RooExtendPdf *sigPdfE = new RooExtendPdf ("sigPdfE", "extended signal 2D PDF", *sigPdf, nSig);
  RooExtendPdf *bkgPdfE = new RooExtendPdf ("bkgPdfE", "extended background 2D PDF", *bkgPdf, nBkg);

  RooAddPdf *fullPdf = new RooAddPdf ("fullPdf", "full 2D PDF", RooArgList (*sigPdfE, *bkgPdfE), RooArgList (nSig, nBkg));
  //  RooProdPdf *fullPdf = new RooProdPdf ("fullPdf", "full 2D PDF", RooArgList (*SigmaPdf, *K0SPdf));

  RooFitResult *fitRes = fullPdf->fitTo (dh, SumW2Error(kTRUE),  Minimizer("Minuit"), PrintLevel(1),
                      Save());

  RooPlot* frame_ppiz = mppi0.frame(Title("Imported TH1 with Poisson error bars")) ;
  RooPlot* frame_pizpiz = mpi0pi0.frame(Title("Imported TH1 with Poisson error bars")) ;
  dh.plotOn(frame_ppiz) ;
  dh.plotOn(frame_pizpiz) ;

  //SigmaPdf->plotOn(frame_ppiz);
  //K0SPdf  ->plotOn(frame_pizpiz);

  fullPdf->plotOn (frame_ppiz);
  fullPdf->plotOn (frame_ppiz, Components(*bkgPdf), LineStyle(kDashed));

  fullPdf->plotOn (frame_pizpiz);
  fullPdf->plotOn (frame_pizpiz, Components(*bkgPdf), LineStyle(kDashed));

  const string cname = std_ext::formatter() << "c" << globalBin;
  TCanvas *c1 = new TCanvas(cname.c_str(), cname.c_str(), 600, 600);
  globalBin++;
  c1->Divide(2,2);

  c1->cd(1);
  frame_ppiz->Draw();

  c1->cd(2);
  frame_pizpiz->Draw();


  // Create and fill ROOT 2D histogram (8x8x8 bins) with contents of dataset
  TH1* hh_data3 = dh.createHistogram("hh_data3",
                     mppi0,Binning(24),
                     YVar(mpi0pi0, Binning(21))) ;


  // Create and fill ROOT 2D histogram (20x20x20 bins) with sampling of pdf
  TH1* hh_pdf3 = fullPdf->createHistogram("hh_model3",mppi0,Binning(100),YVar(mpi0pi0,Binning(100))) ;
    hh_pdf3->SetFillColor(kBlue) ;

  c1->cd(3);
  hh_data3->Draw("COLZ");


  c1->cd(4);
  hh_pdf3->Draw("COLZ");

//    cout << "*** Found nSig = " << nSig.getVal() << " ± " << nSig.getError() << endl;
  if (fitRes->status() == 0)
      return {nSig.getVal(),nSig.getError()};
  return {NaN,NaN};
}

vector<TH2D*> getHists(WrapTFileInput& f, const string& hpath,const string& sndpath, const string& histname, const unsigned nbins)
{
    vector<TH2D*> hs;
    for (auto egbin = 0u ; egbin < nbins; ++egbin)
    {
        TH2D* h = nullptr;

        if(!f.GetObject(std_ext::formatter() << hpath << egbin << "/" << sndpath << histname, h)) {
            LOG(FATAL) << "Cannot find " << hpath;
            continue;
        };
        hs.push_back(h);
    }
    return hs;
};

size_t FillGraphErrors(TGraphErrors* graph, const double x, const double y, const double xerr, const double yerr)
{
    auto N = graph->GetN();
    graph->SetPoint(N,x,y);
    graph->SetPointError(N,xerr,yerr);
    return graph->GetN();
}



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

    auto cmd_nEgBins       = cmd.add<TCLAP::ValueArg<unsigned>>("","nEgBins","Number of bins in Egamma (should match final plot)",false,10,"unsigned");


    auto cmd_output        = cmd.add<TCLAP::ValueArg<string>>("o","output","Output file",false,"","filename");


    TCLAP::ValuesConstraintExtra<decltype(ExpConfig::Setup::GetNames())> allowedsetupnames(ExpConfig::Setup::GetNames());
    auto cmd_setup = cmd.add<TCLAP::ValueArg<string>>("s","setup","Setup to determine tagged photon energy bins",false,"Setup_2014_10_EPT_Prod",&allowedsetupnames);



    cmd.parse(argc, argv);
    if(cmd_verbose->isSet()) {
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());
    }


    ExpConfig::Setup::SetByName(cmd_setup->getValue());

    argc=0;
    TRint app("SigmaPlus_fit",&argc,argv,nullptr,0,true);

    unique_ptr<WrapTFileOutput> outFile;
    if(cmd_output->isSet()) {
        // cd into masterFile upon creation
        outFile = std_ext::make_unique<WrapTFileOutput>(cmd_output->getValue(), true);
    }
    HistogramFactory histfac("SigmaPlus_fit");


    const auto nEgammaBins = cmd_nEgBins->getValue();
    auto tagger = ExpConfig::Setup::GetDetector<TaggerDetector_t>();
    const auto eMin = tagger->GetPhotonEnergy(tagger->GetNChannels()-1);
    const auto eMax = tagger->GetPhotonEnergy(0);
    const double binwidth = (eMax - eMin )/ nEgammaBins;

    WrapTFileInput input_plotter(cmd_PlotFile->getValue());
    auto dalitzHists = getHists(input_plotter,cmd_HistPath->getValue(), cmd_relpath_data->getValue(), cmd_HistName->getValue(), nEgammaBins);
    WrapTFileInput input_plotter_mc(cmd_MCTrue->getValue());
    auto dalitzHists_MC = getHists(input_plotter_mc, cmd_HistPath->getValue(), cmd_relpath_mc->getValue(), cmd_HistName->getValue(), nEgammaBins);

    auto dataGraph = histfac.makeGraphErrors("counts data","nData");
    auto mcGraph = histfac.makeGraphErrors("counts mc","nMC");
    auto finalGraph = histfac.makeGraphErrors("cross section","finalPlot");

    for (auto i = 0u ; i < dalitzHists.size() ; ++i)
    {
        const auto result = fitHist(dalitzHists.at(i));
        const auto resultmc = fitHist(dalitzHists_MC.at(i));
        const auto resultNorm = result / resultmc;
        if (!isnan(result.v) && !isnan(resultmc.v) )
        {
            const auto e = eMin + binwidth * ( .5 + i);
            FillGraphErrors(dataGraph, e,result.v, 0 ,result.e);
            FillGraphErrors(mcGraph, e,resultmc.v, 0 ,resultmc.e);
            FillGraphErrors(finalGraph, e,resultNorm.v, 0 ,resultNorm.e);
        }
    }

    auto cdata = new TCanvas("result data","result data",600,600);
    dataGraph->Draw("AP");
    auto cmc = new TCanvas("result mc","result mc",600,600);
    mcGraph->Draw("AP");
    auto cfinal = new TCanvas("result","result",600,600);
    finalGraph->Draw("AP");

    auto setLabels = [](TGraphErrors* g)
    {
        g->GetXaxis()->SetTitle("E_{#gamma} [MeV]");
        g->GetYaxis()->SetTitle("# events");
    };
    for (auto g: {dataGraph,mcGraph,finalGraph})
        setLabels(g);



    if(outFile)
        LOG(INFO) << "Close ROOT properly to write data to disk.";

    app.Run(kTRUE); // really important to return...
    if(outFile)
        LOG(INFO) << "Writing output file...";
    outFile = nullptr;

    return 0;
}
