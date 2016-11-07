#include "TwoPi0_MCSmearing_tools.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "root-addons/analysis_codes/Math.h"

#include <iostream>

using namespace std;
using namespace ant;

PeakFitResult_t ant::TowPi0_MCSmearing_Tool::Fit(TH1* h)
{

    const double r_min = 80.0;
    const double r_max = 200.0;
    const double exp_center = 135.0;
    const double exp_width= 12.0;
    const int npx= 300;

    TF1* sig = new TF1("sig", "gaus", r_min, r_max);
    sig->SetLineColor(kGreen);
    sig->SetNpx(npx);

    // height
    sig->SetParameter(0, 0.5 * h->GetMaximum());

    sig->SetParameter(1, exp_center);

    // width
    sig->SetParameter(2, 12.0);



    TF1* bg = new TF1("bg", "pol3", r_min, r_max);
    bg->SetLineColor(kBlue);

    bg->SetParameter(0,0);
    bg->SetParName(0, "BG p_{0}");
    bg->SetParameter(1,0);
    bg->SetParName(1, "BG p_{1}");
    bg->SetParameter(2,0);
    bg->SetParName(2, "BG p_{2}");
    bg->SetParameter(3,0);
    bg->SetParName(3, "BG p_{3}");

    TFSum::FitRanged(h, bg, r_min, exp_center - 2*exp_width, exp_center+2*exp_width, r_max);
    //bg->FixParameter(0, bg->GetParameter(0));
    //bg->FixParameter(1, bg->GetParameter(1));
    //bg->FixParameter(2, bg->GetParameter(2));


    TFSum* sum = new TFSum("sum", sig, bg, r_min, r_max);
    sum->SetNpx(npx);


    TCanvas* c = new TCanvas();
    c->SetTitle(Form("Fit to %s", h->GetTitle()));
    h->SetStats(true);
    gStyle->SetOptFit(1);
    h->Draw();
    h->Fit(sum->Function(), "REM0NB");

    sum->SyncToFcts();

    sum->Draw();


    const double total_area = sum->Function()->Integral(r_min, r_max);
    const double bg_area    =  bg->Integral(r_min, r_max);
    const double sig_area   = total_area - bg_area;

    const double sig_to_bg = sig_area / bg_area;
    const double peak_pos  = sig->GetParameter(1);

    cout << "Mass offset = " << peak_pos - exp_center << " MeV\n";
    cout << "Sig/BG      = " << sig_to_bg << "\n";
    cout << "Sig         = " << sig_area << endl;

    //TODO return chi2/dof
    return PeakFitResult_t(0, peak_pos, sig->GetParameter(2));

}
