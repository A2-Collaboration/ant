#include "Fits.h"

#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "root-addons/analysis_codes/Math.h"
#include "base/ParticleType.h"
#include "TStyle.h"
#include "TLatex.h"
#include "base/std_ext/string.h"
#include "TGraph.h"

using namespace ant;
using namespace std;

Fits::FitResult Fits::FitEtaCalib(TH1* h, const double r_min, const double r_max) {
    return FitPeakPol4(h, ParticleTypeDatabase::Eta.Mass(), 21.0, r_min, r_max);
}

Fits::FitResult Fits::FitPi0Calib(TH1* h, const double r_min, const double r_max) {
    return FitPeakPol4(h, ParticleTypeDatabase::Pi0.Mass(), 11.0, r_min, r_max);
}

Fits::FitResult Fits::FitPeakPol4(TH1* h, const double mass, const double expected_width, const double r_min, const double r_max) {

    const int    npx   = 500;

    TF1* sig = new TF1("sig", "gaus", r_min, r_max);
    sig->SetLineColor(kGreen);
    sig->SetNpx(npx);

    // height
    sig->SetParameter(0, 0.5 * h->GetMaximum());

    // peak position
    sig->SetParameter(1, mass);

    // width
    sig->SetParameter(2, expected_width);



    TF1* bg = new TF1("bg", "pol4", r_min, r_max);
    bg->SetLineColor(kBlue);

    bg->SetParameter(0,0);
    bg->SetParName(0, "BG p_{0}");
    bg->SetParameter(1,0);
    bg->SetParName(1, "BG p_{1}");
    bg->SetParameter(2,0);
    bg->SetParName(2, "BG p_{2}");
    bg->SetParameter(3,0);
    bg->SetParName(3, "BG p_{3}");
    bg->SetParameter(4,0);
    bg->SetParName(4, "BG p_{4}");

    const auto peak_range = interval<double>::CenterWidth(mass,2*2*mass);

    TFSum::FitRanged(h, bg, r_min, peak_range.Start(), peak_range.Stop(), r_max);


    TFSum* sum = new TFSum("sum", sig, bg, r_min, r_max);
    sum->SetNpx(npx);


    //TCanvas* c = new TCanvas();
    //c->SetTitle(Form("Fit to %s", h->GetName()));
    h->SetStats(true);
    gStyle->SetOptFit(1);
    //h->Draw();
    h->Fit(sum->Function(), "REM0NB");

    sum->SyncToFcts();

    sum->Draw();


    const double total_area = sum->Function()->Integral(r_min, r_max);
    const double bg_area    =  bg->Integral(r_min, r_max);
    const double sig_area   = total_area - bg_area;
    const double peak_pos  = sig->GetParameter(1);
    const double peak_width = sig->GetParameter(2);
    const double relwidth   = peak_width / peak_pos;

    cout << "Mass offset = " << peak_pos - mass << " MeV\n";
    cout << "sigma/pos   = " <<  relwidth << "\n";
    cout << "Sig         = " << sig_area << endl;

    const string text = std_ext::formatter()
            << "N = " << sig_area << " "
            << "#sigma/#mu = " << relwidth;
    auto l = new TLatex(peak_pos, sum->Function()->Eval(peak_pos), text.c_str());
    l->Draw();


    return FitResult(peak_pos, sig_area, sig->GetParameter(2), 0);

}

void Fits::FitSlicesPi0(TH2 *h2)
{
    double minEnergy=150;
    double maxEnergy=450;
    TGraph* g1 = new TGraph();
    TGraph* g1_rel = new TGraph();
    int k=0;
    for(int i=1; i>0; ++i) {
        TH1* b = h2->ProjectionX(Form("x%d",i),i,i+1);
        double e = h2->GetYaxis()->GetBinCenter(i);

        if (e < maxEnergy && e > minEnergy)
        {
            new TCanvas();
            b->Draw();
            auto result = FitPi0Calib(b);
            g1->SetPoint(k,e,result.pos);
            g1_rel ->SetPoint(k,e,(result.pos/135-1) * 100);
            k++;
        }
        else {
            k = 0;
        }
        if(e > maxEnergy){
            break;
        }
    }
    new TCanvas();
    g1->Draw();
    g1->SetTitle("Position of the pi0 peak in different energy intervals of 25 MeV");
    g1->GetXaxis()->SetTitle("Energy of the photons [MeV]");
    g1->GetYaxis()->SetTitle("Position of pi0 peak [MeV]");

    new TCanvas();
    g1_rel->Draw();
    g1_rel->SetTitle("Position of the pi0 peak in different energy intervals of 25 MeV");
    g1_rel->GetXaxis()->SetTitle("Energy of the photons [MeV]");
    g1_rel->GetYaxis()->SetTitle("Deviation from the 135 MeV peak [%] ");
}


