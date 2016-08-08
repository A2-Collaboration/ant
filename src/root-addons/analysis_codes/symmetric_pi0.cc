#include "symmetric_pi0.h"

#include "TTree.h"
#include "TH2.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TCut.h"

#include "TreeTools.h"

using namespace std;
using namespace ant;

void SymmetricPi0::Analyse(TTree* tree)
{
    new TCanvas();

    const TCut cb = "Theta1*TMath::RadToDeg()>22 && Theta2*TMath::RadToDeg()>22";
 //   const TCut ERange = "E1/(E1-E)<.01";

    tree->Draw("IM:E>>a(1000,0,1000,120,80,200)", cb , "colz");

    TH2* a = NULL;
    gDirectory->GetObject("a",a);
    if(!a) {
        cout << "a not found!" << endl;
        return;
    }

    findMax(a);
}

void SymmetricPi0::findMax(TH2* h) {

    const int start = 80;
    const int end   = 500;
    const int width = 50;
    const int n = (end-start)/width +1;

    TGraphErrors* g= new TGraphErrors(n);
    TGraphErrors* sigmas= new TGraphErrors(n);

    TH1D* pr = NULL;

    TF1* f = new TF1("f", "gaus+pol3(3)", 80, 200);
    f->SetParameter(0, 6.70979e+03);
    f->SetParameter(1, 1.34825e+02);
    f->SetParLimits(1, 120, 150);
    f->SetParameter(2, 1.12354e+01);
    f->SetParLimits(2, 1,50);
    f->SetParameter(3, 4.55534e+02);
    f->SetParameter(4, 7.29161e+01);
    f->SetParameter(6, -3.24012e-01);
    f->SetParameter(7, -1.95152e-05);

    for(int i=0; i<n; ++i) {
        pr = h->ProjectionY(Form("pr_%d",i), start+i*width, start+(i+1)*width);
        pr->SetTitle(Form("Energy Range: %d to %d MeV", start+i*width, start+(i+1)*width));
        pr->SetXTitle("IM [MeV]");
        pr->Draw();
        pr->Fit(f, "R");
        g->SetPoint(i, start+i*width, f->GetParameter(1));
        g->SetPointError(i, 0, f->GetParError(1));
        sigmas->SetPoint(i, start+i*width, f->GetParameter(2));
        sigmas->SetPointError(i, 0, f->GetParError(2));
        gPad->SaveAs(Form("Symmetric_pi0_fit_%d.png",i));
        gPad->SaveAs(Form("Symmetric_pi0_fit_%d.pdf",i));
    }

    new TCanvas();
    g->Draw();

    new TCanvas();
    sigmas->Draw();
}
