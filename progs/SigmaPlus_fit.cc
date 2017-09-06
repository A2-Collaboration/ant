#include "base/Logger.h"

#include "tclap/CmdLine.h"
#include "tclap/ValuesConstraintExtra.h"
#include "base/WrapTFile.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TFile.h"

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



void fitHist (TH2D* histDalitz) {

  RooRealVar mppiz("mppiz","m(p #pi^{0})",1150, 1230) ;
  RooRealVar mpizpiz("mpizpiz","m(#pi^{0} #pi^{0})", 420, 560) ;

  RooDataHist dh ("dh", "dh", RooArgSet(mppiz, mpizpiz), histDalitz);


  RooRealVar mean_Sigma("mean_Sigma","m_{#Sigma}",1188, 1170, 1220);
  RooRealVar sigma_Sigma("sigma_Sigma","#sigma_{#Sigma}",10, 0, 50);
  RooRealVar alpha_Sigma("alpha_Sigma","#alpha_{#Sigma}",2);
  RooRealVar n_Sigma("n_Sigma","n_Sigma",1);
  RooCBShape *theSigma = new RooCBShape("theSigma","crystal ball PDF",mppiz,mean_Sigma,sigma_Sigma,alpha_Sigma,n_Sigma);

  RooRealVar p0_Sigma ("p0_Sigma", "p0_Sigma", 0, -1, 1);
  RooRealVar p1_Sigma ("p1_Sigma", "p1_Sigma", 0, -1, 1);
  RooRealVar p2_Sigma ("p2_Sigma", "p2_Sigma", 0, -1, 1);
  RooRealVar p3_Sigma ("p3_Sigma", "p3_Sigma", 0);

  RooRealVar frac_peak ("frac_peak", "frac_peak", 0.5, 0, 1);

  RooChebychev *bkg_Sigma = new RooChebychev ("bkg_Sigma", "bkg_Sigma",
                          mppiz, RooArgList(p0_Sigma, p1_Sigma, p2_Sigma, p3_Sigma));

  //  RooAddPdf *SigmaPdf = new RooAddPdf ("SigmaPdf", "SigmaPdf",
  //				       RooArgList (*theSigma, *bkg_Sigma),
  //				       RooArgList (frac_peak));

  RooRealVar mean_K0S("mean_K0S","m_{K_{S}}",490, 460, 500);
  RooRealVar sigma_K0S("sigma_K0S","#sigma_{K_{S}}",12, 0, 50);
  RooRealVar alpha_K0S("alpha_K0S","#alpha_{K_{S}}",0.97);
  RooRealVar n_K0S("n_K0S","n_K0S",10);
  RooCBShape *theK0S = new RooCBShape("theK0S","crystal ball PDF for K0S",mpizpiz,mean_K0S,sigma_K0S,alpha_K0S,n_K0S);

  RooRealVar p0_K0S ("p0_K0S", "p0_K0S", 0, -1, 1);
  RooRealVar p1_K0S ("p1_K0S", "p1_K0S", 0, -1, 1);
  RooRealVar p2_K0S ("p2_K0S", "p2_K0S", 0, -1, 1);

  RooRealVar frac_peak_K0S ("frac_peak_K0S", "frac_peak_K0S", 0.5, 0, 1);

  RooChebychev *bkg_K0S = new RooChebychev ("bkg_K0S", "bkg_K0S",
                        mpizpiz, RooArgList(p0_K0S, p1_K0S, p2_K0S));


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

  RooPlot* frame_ppiz = mppiz.frame(Title("Imported TH1 with Poisson error bars")) ;
  RooPlot* frame_pizpiz = mpizpiz.frame(Title("Imported TH1 with Poisson error bars")) ;
  dh.plotOn(frame_ppiz) ;
  dh.plotOn(frame_pizpiz) ;

  //SigmaPdf->plotOn(frame_ppiz);
  //K0SPdf  ->plotOn(frame_pizpiz);

  fullPdf->plotOn (frame_ppiz);
  fullPdf->plotOn (frame_ppiz, Components(*bkgPdf), LineStyle(kDashed));

  fullPdf->plotOn (frame_pizpiz);
  fullPdf->plotOn (frame_pizpiz, Components(*bkgPdf), LineStyle(kDashed));

  TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
  c1->Divide(2,2);

  c1->cd(1);
  frame_ppiz->Draw();

  c1->cd(2);
  frame_pizpiz->Draw();


  // Create and fill ROOT 2D histogram (8x8x8 bins) with contents of dataset
  TH1* hh_data3 = dh.createHistogram("hh_data3",
                     mppiz,Binning(24),
                     YVar(mpizpiz, Binning(21))) ;


  // Create and fill ROOT 2D histogram (20x20x20 bins) with sampling of pdf
  TH1* hh_pdf3 = fullPdf->createHistogram("hh_model3",mppiz,Binning(100),YVar(mpizpiz,Binning(100))) ;
  //  hh_pdf3->SetFillColor(kBlue) ;

  c1->cd(3);
  hh_data3->Draw("COLZ");


  c1->cd(4);
  hh_pdf3->Draw("COLZ");

  if (fitRes->status() == 0)
    cout << "*** Found nSig = " << nSig.getVal() << " ± " << nSig.getError() << endl;
  else
    cout << "*** Fit Failed.  Results unreliable " << endl;


}

template <typename T>
T* getHist(WrapTFileInput& f, const string& hpath) {
    T* h;
    if(!f.GetObject(hpath, h)) {
        LOG(FATAL) << "Cannot find " << hpath;
    };
    return h;
};

int main(int argc, char** argv) {
    SetupLogger();

    TCLAP::CmdLine cmd("SinglePi0_fit", ' ', "0.1");

    auto cmd_verbose       = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"int");

    auto cmd_PlotFile          = cmd.add<TCLAP::ValueArg<string>>("","plotFile","Data input from  singlePi0-Plot",true,"","rootfile");

    auto cmd_HistPath          = cmd.add<TCLAP::ValueArg<string>>("","histName","Data input from  singlePi0-Plot",false,"ppi0_2pi0","rootfile");


    cmd.parse(argc, argv);
    if(cmd_verbose->isSet()) {
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());
    }


    argc=0;
    TRint app("SigmaPlus_fit",&argc,argv,nullptr,0,true);

    WrapTFileInput input_lumi(cmd_PlotFile->getValue());
    auto histDalitz = getHist<TH2D>(input_lumi,cmd_HistPath->getValue());

    fitHist(histDalitz);

    app.Run(kTRUE);

    return 0;
}
