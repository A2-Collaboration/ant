#include "analysis/plot/RootDraw.h"
#include "base/Detector_t.h"
#include "base/Logger.h"
#include "base/PlotExt.h"
#include "base/std_ext/math.h"
#include "base/std_ext/string.h"
#include "base/std_ext/system.h"
#include "base/std_ext/time.h"
#include "base/WrapTFile.h"
#include "expconfig/setups/Setup.h"
#include "expconfig/ExpConfig.h"
#include "expconfig/detectors/Tagger.h"
#include "calibration/DataManager.h"
#include "calibration/modules/TaggEff.h"
#include "tree/TCalibrationData.h"
#include "detail/taggEffClasses.h"
#include "tclap/CmdLine.h"
#include "tree/TAntHeader.h"
#include "calibration/fitfunctions/FitGausPol0.h"
#include "calibration/fitfunctions/FitFunction.h"
#include "calibration/fitfunctions/BaseFunctions.h"
#include "analysis/plot/HistogramFactory.h"

#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TRint.h"
#include "TTree.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include <TMath.h>
#include <TMinuit.h>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"

using namespace ant;
using namespace std;
using namespace ant::analysis;
using namespace ant::std_ext;
using namespace ant::calibration;

// some other bools used in various methods
static volatile bool interrupt = false;
static bool noStore = false;
static bool histOut = false;

// global histograms
HistogramFactory *histfac;
TH1D *hDataEnh = NULL;
TH1D *hDataRawEnh = NULL;
TH1D *hCalcEnh = NULL;
TH1D *hPol = NULL;
//these are really just convenient arrays - they don't ever get plotted.
TH1F *weightHist  = NULL;
TH2F *thetaWeight = NULL;
TH2F *thetaPol    = NULL;
TH2F *thetaTtot   = NULL;
TH2F *thetaItot   = NULL;

// functions
Double_t efit(const Double_t *);
void enhFromHuman(Double_t beamMeV = 883.0, Double_t edgeMeV = 145.0, Double_t spreadMeV = 10.0, Double_t colliDist_m = 2.5, Double_t colliRad_mm = 1.0, Int_t nVec = 4);
void parFromHuman(Double_t beamMeV = 883.0, Double_t edgeMeV = 145.0, Double_t spreadMeV = 10.0, Double_t colliDist_m = 2.5, Double_t colliRad_mm = 1.0, Int_t nVec = 4, Double_t *par=NULL);
void enhFromParams(Double_t *par=NULL);
void init();


//Some enumerators and names
enum {
  THETA,  // [0] theta      main angle responsible for coherent edge cutoffs
  SIGMA,  // [1] sigma      smearing of theta
  THETAR, // [2] thetar     relative angle resonsible for colli cutoffs
  SIGMAR, // [3] sigmar     smearing of colli cutoff angle
  E0MEV,  // [4] beam energy
  NVEC,   // [5] nvec       no of vectors contributing
  IVEC};  // [6] ivec[]     array of intensities of vectors up to nvec.


// Some basic consts etc first
// Consts are all UPPER CASE
//Approx Form factor is F(g^2) = (q^2 + b^(-2)) ^ -2
//Where b= 111 x Z^(-1/3) (x 925 to get into units of crystal lattice)
const Double_t B = 0.247892436;  //where did I get that ? Timm ?
const Double_t A=0.03;           //made up for now, need to get the actual no for this later
const Double_t k=26.5601;        //put in formula for k later (my own stonehenge paper)

const Int_t VECTORS[]={2,4,6,8,10};    //list of the vectors to be included (022,044);

//THESE NEED TO BE CHANGED FOR EACH SETTING (ie comment in/out)
Int_t THETASTEPS = 201;          //no of steps in convoluting with gaussian
Double_t LOWFIT = 20.0;         //how far below the peak in MeV to start fitting

//-- Tagger energy range
vector<double> energy_bin_edges;
int firsttch = 0;

//-- Some stuff needed for the fit
Double_t fitMinEnergy;
Int_t fitMinBin;
Double_t fitMaxEnergy;
Int_t fitMaxBin;
Int_t verbose=0;
Double_t bestPar[10];
Double_t bestChisq;
TF1 *gausFit;
Bool_t isInit=kFALSE;

auto failExit = [] (const string& message)
{
    LOG(ERROR) << message;
    exit(EXIT_FAILURE);
};

int main( int argc, char** argv )
{
    SetupLogger();

    signal(SIGINT, [] (int) {
        LOG(INFO) << ">>> Interrupted";
        interrupt = true;
    });

    //-- Read from the command line
    TCLAP::CmdLine cmd("Ant-makeLinPol - Create LinPol tables from tagging efficiency runs - one polarised and one unpolarised run", ' ', "0.1");
    //--- other settings:
    auto cmd_output    = cmd.add<TCLAP::ValueArg<string>>("o", "output", "Output file", false, "", "filename");
    auto cmd_colliDist = cmd.add<TCLAP::ValueArg<double>>("", "colliDist", "Distance to collimeter? (m)",false,3.0,"colliDist");
    auto cmd_colliRad  = cmd.add<TCLAP::ValueArg<double>>("", "colliRad", "Radius of collimeter? (mm)",false,1.0,"colliRad");
    auto cmd_nVec      = cmd.add<TCLAP::ValueArg<int>>("", "nVec", "Number of harmonics????",false,3,"nVec");
    //--- switches
    auto cmd_batchmode = cmd.add<TCLAP::SwitchArg>("b", "batch",   "Run in batch mode (no ROOT shell afterwards)");
    auto cmd_nostore   = cmd.add<TCLAP::SwitchArg>("n", "nostore", "Don't store polarisation tables in the calibration database, only show results");
    //--- files
    auto cmd_filelist  = cmd.add<TCLAP::UnlabeledMultiArg<string>>("inputfiles", "Input files to read from", true, "Pol-inputfile UnPol-inputfile");
    cmd.parse(argc, argv);
    auto fileList = cmd_filelist->getValue();
    auto colliDist_m = cmd_colliDist->getValue();
    auto colliRad_mm = cmd_colliRad->getValue();
    auto nVec = cmd_nVec->getValue();
    noStore = cmd_nostore->isSet();
    histOut = cmd_output->isSet();

    //-- what the heck does this do?
    argc=0; // prevent TRint to parse any cmdline (?)
    TRint app("Ant-makeLinPol",&argc,argv,nullptr,0,true);
    unique_ptr<WrapTFileOutput> masterFile;
    if(cmd_output->isSet()) {
        // cd into masterFile upon creation
        masterFile = std_ext::make_unique<WrapTFileOutput>(cmd_output->getValue(), true);
    }

    //-- Open the inputfiles, fetch the scalar distributions and create the enhancement spectrum
    TFile *fPol = new TFile((TString)fileList.at(0),"read");
    TFile *fUnpol = new TFile((TString)fileList.at(1),"read");
    TH1D *hPolScaler_in = (TH1D*)fPol->Get("ProcessTaggEff/scalerHits");
    TH1D *hUnpolScaler_in = (TH1D*)fUnpol->Get("ProcessTaggEff/scalerHits");
    //--- Get the setup name from the first root-file
    TAntHeader* header; fPol->GetObject("AntHeader",header);
    ExpConfig::Setup::SetByName(header->SetupName);

    //-- Get some more info from the setup
    auto Tagger = ExpConfig::Setup::GetDetector<expconfig::detector::Tagger>();
    if (!Tagger) throw std::runtime_error("No Tagger found");
    unsigned nchannels = Tagger->GetNChannels();
    double beamMeV = ExpConfig::Setup::Get().GetElectronBeamEnergy();

    //-- Create histograms
    auto bins_tagger = BinSettings(nchannels);
    histfac = new HistogramFactory("makeLinPol");
    TH1D *hPolScaler = histfac->makeTH1D("PolScaler","Tagger channel","",bins_tagger,"hPolScaler",true);
    TH1D *hUnpolScaler = histfac->makeTH1D("UnpolScaler","Tagger channel","",bins_tagger,"hUnpolScaler",true);
    //---- Create an x-axis binsetting of the photon energy
    //----- Start the energy range from the first tagger channel used and stop at the last
    //      Currently assuming that one is only turning off channels in the low region and using all up to highest!!! - should be generalised!
    for(int t=0; t<(int)nchannels; t++){
        if(Tagger->HasElementFlags(t, Detector_t::ElementFlag_t::Broken)) continue;
        else {firsttch=t; break;}
    }
    unsigned nchannelsused = nchannels - firsttch;
    vector<double> photon_energies;
    for(unsigned ch=firsttch;ch<nchannels;ch++)
        photon_energies.push_back(Tagger->GetPhotonEnergy(ch));
    //----- make some bin edges half way between the energy values.
    vector<double> energy_bins(nchannelsused+1);
    for(unsigned b=0;b<(nchannelsused-1);b++){
      energy_bins.at(b+1)=0.5*(photon_energies.at(b)+photon_energies.at(b+1));
    }
    //----- and the top and bottom have width of the adjacent bin
    energy_bins.at(0) = energy_bins.at(1) - (energy_bins.at(2) - energy_bins.at(1));
    energy_bins.at(nchannelsused) = energy_bins.at(nchannelsused-1) + (energy_bins.at(nchannelsused-1)- energy_bins.at(nchannelsused-2));
    for(unsigned i=0; i<nchannelsused+1;i++){
        energy_bin_edges.push_back(energy_bins[i]);
    }
    hDataEnh = histfac->makeTH1D("DataEnhancement",energy_bins,"Photon energy","","hDataEnh",true);
    hDataRawEnh = histfac->makeTH1D("DataRawEnhancement",energy_bins,"Photon energy","","hDataRawEnh",true);

    //--- And fill the data histograms
    hPolScaler->Add(hPolScaler_in);
    hUnpolScaler->Add(hUnpolScaler_in);
    double polbc = 0, unpolbc=0, enh=0, enherr=0;
    for(unsigned i=1; i<nchannelsused+1;i++){
        polbc = hPolScaler->GetBinContent(i+firsttch);
        unpolbc = hUnpolScaler->GetBinContent(i+firsttch);
        if(polbc == 0 || unpolbc == 0){
            enh = 0; enherr = 0;
        }
        else {
            enh = polbc/unpolbc;
            enherr = sqrt(polbc/pow(unpolbc,2) + pow(polbc,2)/pow(unpolbc,3));
        }
        hDataEnh->SetBinContent(i,enh);
        hDataRawEnh->SetBinContent(i,enh);
        hDataEnh->SetBinError(i,enherr);
        hDataRawEnh->SetBinError(i,enherr);
    }

    //-- Some stuff needed later on
    Double_t par[10];

    //-- In cbremFit_R there would here be some loops to get rid of zeros and spikes (aka THE SMOOTHY BIT)

    //-- Find a reasonable minumum spot to set to 1 for the baseline.
    double lowmean = 1000000.0;
    for(int n=(firsttch+(int)(0.05*(float)nchannels));n<=(int)(0.75*(float)nchannels);n++){
        if((hDataEnh->GetBinContent(n)>0.0)&&(hDataEnh->GetBinContent(n-2)>0.0)&&(hDataEnh->GetBinContent(n+2)>0.0)&&(hDataEnh->GetBinContent(n-1)>0.0)&&(hDataEnh->GetBinContent(n+1)>0.0)){
            if((hDataEnh->Integral(n-2,n+2)<lowmean)){
                lowmean=hDataEnh->Integral(n-2,n+2);
            }
        }
    }
    // why this scaling factor?
    hDataEnh->Scale(5./(lowmean));
    cout << "Low Mean = " << lowmean << endl;

    //-- Now try to make some guesses at the initial parameters
    gausFit = functions::GausPol<0>::getTF1();
    gausFit->SetNpx(1000);
    gausFit->SetRange(hDataEnh->GetBinCenter(hDataEnh->GetMaximumBin()),hDataEnh->GetBinCenter(hDataEnh->GetMaximumBin())+100.0);
    gausFit->SetParameter(1,hDataEnh->GetBinCenter(hDataEnh->GetMaximumBin()));
    gausFit->SetParameter(2,10.0);
    gausFit->SetParameter(3,1.0);
    hDataEnh->Fit(gausFit,"rN");
    lowmean=0.0;
    double fitedge=0.0;
    //Get the edge from the derivative
    for(float d = hDataEnh->GetBinCenter(hDataEnh->GetMaximumBin());d < hDataEnh->GetBinCenter(hDataEnh->GetMaximumBin()+90.0);d+=0.1){
        if(gausFit->Derivative(d)<lowmean){
            lowmean=gausFit->Derivative(d);
            fitedge=d;
        }
    }
    cout << "edge = " << fitedge << " MeV" << endl;

    //Now we have enough information to set the basic parameters
    parFromHuman(beamMeV,fitedge,gausFit->GetParameter(2),colliDist_m,colliRad_mm,nVec,par);

    //set the intensities
    for(int v=0;v<par[NVEC];v++){                                               //give the vectors intensities
      par[IVEC+v] = hDataEnh->GetMaximum()*2.0/((Double_t)VECTORS[v]*(Double_t)VECTORS[v]);      //tailing off as 1/VECTORS[v]^2
      //cout << IVEC+v << "  v   " << par[IVEC+v] << endl;
    }

    enhFromParams(par);

    //Redo the intensities according to a the calc / data ration
    double scalefac=hDataEnh->GetMaximum()/hCalcEnh->GetMaximum();
    for(int v=0;v<par[NVEC];v++){                                               //give the vectors intensities
      par[IVEC+v]*=scalefac;
    }
    enhFromParams(par);

    //-- ?? DRAW SOMETHING HERE? ??

    //-- Time to fit

    //--- Set the range of the fit to be some sensible amount below peak and just past the 2nd peak (that calculation is a bit mysterious o_O ).
    fitMinBin = hDataEnh->FindBin(hDataEnh->GetBinCenter(hDataEnh->GetMaximumBin())-LOWFIT);
    fitMaxBin = hDataEnh->FindBin( (par[E0MEV]/(( ( ( 2.0/4.0 )*(( par[E0MEV]/hDataEnh->GetBinCenter(hDataEnh->GetMaximumBin()) ) -1.0 ) ) +1.0 )) ) + 70.0 );
    //--- Create the minuit fitter
    ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Simple");
    //---- set tolerance , etc...
    min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
    min->SetMaxIterations(10000);  // for GSL
    min->SetTolerance(0.001);
    min->SetPrintLevel(1);
    ROOT::Math::Functor ft(&efit,IVEC+nVec);
    min->SetFunction(ft);
    //---- set the variables
    min->SetLimitedVariable(THETA,   "Theta",   par[THETA],      par[THETA]/100.0,  0.95*par[THETA], 1.05*par[THETA]);
    min->SetLimitedVariable(SIGMA,   "Sigma",   2.5*par[SIGMA],  par[SIGMA]/100.0,  par[SIGMA],  5.0*par[SIGMA]);
    min->SetLimitedVariable(THETAR,  "Thetar",  2.0*par[THETAR],     par[THETAR]/100.0, 0.2*par[THETAR], 5.0*par[THETAR]);
    min->SetLimitedVariable(SIGMAR,  "Sigmar",  0.5*par[SIGMAR], par[SIGMAR]/100.0, 0.1*par[SIGMAR], 20.0*par[SIGMAR]);
    min->SetFixedVariable(E0MEV,     "E0MeV",   par[E0MEV]);  //no of vectors
    min->SetFixedVariable(NVEC,      "Nvec",    par[NVEC]);  //no of vectors
    Char_t name[30];
    for(int n=0;n<nVec;n++){
        sprintf(name,"Vec0%d%d", VECTORS[n],VECTORS[n]);
        min->SetVariable(n+IVEC, name, par[n+IVEC], par[n+IVEC]/100.0);
    }
    //--- do the fit
    bestChisq=100000.00;
    min->Minimize();

    //--- use the resulting parameters
    enhFromParams(bestPar);

    TCanvas *c = new TCanvas("c","c",800,800); c->Divide(2,2);
    c->cd(1); hDataEnh->Draw(); gausFit->Draw("same");
    c->cd(2); hDataRawEnh->Draw();
    c->cd(3); hCalcEnh->Draw();
    c->cd(4); hPol->Draw();


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

    return EXIT_SUCCESS;
}

//The main customized fitting function which gets called by MINUIT
Double_t efit(const Double_t *parms){

    Double_t chisq = 1.0;
    Double_t delta;
    Double_t b1,b2;
    Double_t err;
    Double_t *par = (Double_t*)parms;

    hCalcEnh->Reset("ICE"); //reset the histogram
    //  cout << "par[0]= " << par[0] << endl;

    //call the function to make the enhancement and polarization
    enhFromParams(par);

    chisq = 1.0;
    //loop over all the required bins in the histogram to work out a chisq
    for(int n=fitMinBin;n<=fitMaxBin;n++){
        b1=hCalcEnh->GetBinContent(n);
        b2=hDataEnh->GetBinContent(n);
        err=1.0;
        delta=(b1-b2)/err;
        chisq+=(delta*delta);
        //note - not a proper chisq because its an enhancement
    }

    fprintf(stderr,"Chisq: \t%6.2f\t\r",chisq);

    if(chisq<bestChisq){
        bestChisq=chisq;
        for(int n=0;n<10;n++){
            bestPar[n]=par[n];
        }
    }
    return chisq;
}

void parFromHuman(Double_t beamMeV, Double_t edgeMeV, Double_t spreadMeV, Double_t colliDist_m, Double_t colliRad_mm, Int_t nVec, Double_t *par){

  //takes some physical quantities and makes them into parameters, then calls the
  //enhFromParams function.

  //  Double_t par[10];                                                           //array of parameters
  Int_t g = 2;                                                                //variables used in CLAS note
  Double_t E0 = beamMeV;
  Double_t Eg = edgeMeV;


  par[THETA]  = k/(g*E0*E0*((1/Eg)-(1/E0)));                                  //theta from edge and beam energy
  par[SIGMA]  = (par[THETA]-(k/(g*E0*E0*((1/(Eg-spreadMeV))-(1/E0)))))/3.0;   //spread in theta from spread in edge
  par[THETAR] = E0*0.001*5.0*colliRad_mm/colliDist_m;                         //cut from collimator
  par[SIGMAR] = par[THETAR]*par[SIGMA]/par[THETA];                            //smear in above same fractional sigma as above
  par[E0MEV]  = E0;                                                           //beam energy
  par[NVEC]   = (Double_t)nVec;                                                         //no of harmonics

  for(int v=0;v<par[NVEC];v++){                                               //give the vectors intensities
    par[IVEC+v] = 2.0/(Double_t)VECTORS[v];                                   //tailing off as 1/VECTORS[v]
    //cout << IVEC+v << "  v   " << par[IVEC+v] << endl;
  }
}

void enhFromParams(Double_t *par){
  //make an enhancement and corresponding polarization from some the parameters as defined in the CLAS note.
  //this function is can be called stand alone, but will also ba called many times from the fitting function

  Double_t xd[10];
  Double_t xc[10];
  Double_t Q[10];
  Double_t cohContrib;
  Double_t cohTotal;
  Double_t phiTotal;
  Double_t etotal;
  Double_t ptotal;
  Double_t x=0.0;
  Int_t    g=0;
  Double_t weight=0.0;
  Double_t weightSum=0.0;
  Double_t polSum=0.0;
  Double_t phi,chi,cd;
  Double_t amo;
  Int_t jbin=0;

  //loop over sigma
  // for(int p=0;p<10;p++){
  //  cout << p << ": " << par[p] << ", ";
  //}
  //cout << endl;

  // if needed, make some hists
  int nchannelsused = energy_bin_edges.size() - 1;
  double energy_bins[nchannelsused+1];
  for(int i=0; i<nchannelsused+1;i++){
      energy_bins[i] = energy_bin_edges.at(i);
  }
  if(!hCalcEnh){
      hCalcEnh = histfac->makeTH1D("CalculatedEnhancement",energy_bin_edges,"Photon energy","","hCalcEnh",true);
      hPol = histfac->makeTH1D("Polarisation",energy_bin_edges,"Photon energy","","hPol",true);
      hCalcEnh->SetMinimum(0);
      hPol->SetMinimum(0);
      hPol->SetMaximum(1);
  }
  //-- make some useful arrays/histograms not plotted
  if(!thetaPol){
      weightHist   = new TH1F("weightHist",  "weightHist", THETASTEPS+1, 0, THETASTEPS+1 );
      thetaWeight  = new TH2F("thetaWeight", "thetaWeight",nchannelsused,energy_bins, THETASTEPS+1,0, THETASTEPS+1);
      thetaPol     = new TH2F("thetaPol",    "thetaPol",   nchannelsused,energy_bins, THETASTEPS+1,0, THETASTEPS+1);
      thetaItot    = new TH2F("thetaItot",   "thetaItot",  nchannelsused,energy_bins, THETASTEPS+1,0, THETASTEPS+1);
  }

  //reset them all for fresh filling
  hCalcEnh->Reset("ICE");
  hPol->Reset("ICE");
  thetaPol->Reset("ICE");
  thetaItot->Reset("ICE");
  weightHist->Reset("ICE");
  thetaWeight->Reset("ICE");

  for(Double_t j=par[THETA]-3.0*par[SIGMA];j<=par[THETA]+3.001*par[SIGMA];j+=(6.0*par[SIGMA])/THETASTEPS){
      weight=TMath::Gaus(j,par[THETA],par[SIGMA]);   //get the weight from the gaussian
      weightSum+=weight;                             //add to sum
      //find the discontinuity for each vector
      for(int v=0;v<par[NVEC];v++){
          g=VECTORS[v];
          xd[v]=1.0/((k/(g*par[E0MEV]*j))+1.0);
          Q[v]=(1.0-xd[v])/xd[v];
          xc[v]=xd[v]/(1+((par[THETAR]*par[THETAR])*(1-xd[v])));
      }

      //loop over all bins in the histogram
      for(int bin=1;bin<=hCalcEnh->GetNbinsX();bin++){
          x=hCalcEnh->GetBinCenter(bin)/par[E0MEV];            //find the value of the bin
          amo=1/x;                                    //assume amo = inc = 1/x over regio of interest
          cohTotal=0.0;
          phiTotal=0.0;

          //loop over all the vectors
          for(int v=0;v<par[NVEC];v++){
              if(x>xd[v]) continue;           //only do up to x_dg

              //work out chi and phi
              phi=(2*Q[v]*Q[v]*x*x)/((1-x)*(1+((1-x)*(1-x))-((4*Q[v]*Q[v]*x*x/(1-x))*(((1-x)/(Q[v]*x))-1))));
              chi=((Q[v]*Q[v]*x)/(1-x))*(1+((1-x)*(1-x))-((4*Q[v]*Q[v]*x*x/(1-x))*(((1-x)/(Q[v]*x))-1)));
              //	cout << j  << "  " << chi << endl;
              cd=0.5*(1+TMath::Erf((x-xc[v])/(TMath::Sqrt(2)*par[SIGMAR])));

              //get coherent contrib for the vector
              cohContrib=cd*par[IVEC+v]*chi;

              //add to the total and update the phi total
              cohTotal+=cohContrib;
              phiTotal+=cohContrib*phi;
          }
          //MINEcout<<cohTotal<<" "<<phiTotal<<endl;
          if(cohTotal>0.0) {
              phiTotal/=cohTotal;   //divide by the cohTotal to get the weighted dmean phi
              //cout << x << " " << phiTotal << " " << cohTotal << " " << weight << endl;
          }

          //enhancement = coherent total + inc (or amo).
          etotal=(amo+cohTotal)/amo;
          //and pol like this
          //      ptotal=phiTotal*cohTotal/(cohTotal + amo);
          ptotal=phiTotal*cohTotal;

          //add the weighted contribution to the enhancement
          hCalcEnh->Fill(x*par[E0MEV],weight*etotal);

          //keep the pol for this x,theta coord
          thetaPol->Fill(x*par[E0MEV],jbin,ptotal);

          //keep the total intensity for this x,theta coord
          thetaItot->Fill(x*par[E0MEV],jbin,cohTotal+amo);
      }

      //save the weight for this theta point
      weightHist->Fill(jbin,weight);
      jbin++;
  }
  //normalize the sum of the weighted enhancements
  hCalcEnh->Scale(1.0/weightSum);

  //loop over each x bin, adding the weighted contribs from each theta pos
  for(int bin=1; bin<=hPol->GetNbinsX(); bin++){
      weightSum=0.0;
      polSum=0.0;

      for(int jb=1;jb<=weightHist->GetNbinsX();jb++){
          weight=weightHist->GetBinContent(jb);
          //      polSum+=thetaPol->GetBinContent(bin,jb)*thetaItot->GetBinContent(bin,jb)*weight;
          polSum+=thetaPol->GetBinContent(bin,jb)*weight;
          weightSum+=thetaItot->GetBinContent(bin,jb)*weight;
          //polSum+=thetaPol->GetBinContent(bin,jb)*weight;
          //weightSum+=weight;
      }
      polSum/=weightSum;
      hPol->Fill(hPol->GetBinCenter(bin),polSum);
  }
}

void init(){
    if(isInit) return;
    gausFit = functions::GausPol<0>::getTF1();
    gausFit->SetNpx(1000);
    isInit=kTRUE;
}
