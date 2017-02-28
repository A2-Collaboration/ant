#include "Omega.h"

#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include <iostream>
#include "root-addons/analysis_codes/Math.h"
#include "TAxis.h"
#include "analysis/plot/root_draw.h"
#include "TGraph.h"
#include "base/ParticleType.h"

#include "base/math_functions/CrystalBall.h"

using namespace ant;
using namespace std;

//struct SignalFct_t {

//    TF1* fct = nullptr;

//    virtual double getPos() const =0;
//    virtual double getHeight() const =0;

//    virtual double getArea(const double r_min, const double r_max) const {
//        return fct->Integral(r_min, r_max);
//    }

//    SignalFct_t(TF1* f): fct (f) {};
//};

//struct GausPeak : SignalFct_t {
//    static constexpr double expected_width = 15.0;

//    GausPeak(const double r_min, const double r_max):
//        SignalFct_t(new TF1("sig", "gaus", r_min, r_max)) {
//        fct->SetParameter(0, 0.5 * h->GetMaximum());

//        // position
//        if(fixOmegaMass)
//            fct->FixParameter(1, ParticleTypeDatabase::Omega.Mass());
//        else
//            fct->SetParameter(1, ParticleTypeDatabase::Omega.Mass());

//        // width
//        fct->SetParameter(2, expected_width);
//    }
//};

Omega::FitResult Omega::FitHist(TH1 *h, const bool fixOmegaMass, const double r_min, const double r_max) {

    const int    npx   = 500;
    const double omega_mass     = ParticleTypeDatabase::Omega.Mass();
    const double expected_width =  15.0;

    TF1* sig = new TF1("sig", "gaus", r_min, r_max);
    sig->SetLineColor(kGreen);
    sig->SetNpx(npx);

    // height
    sig->SetParameter(0, 0.5 * h->GetMaximum());

    // position
    if(fixOmegaMass)
        sig->FixParameter(1, omega_mass);
    else
        sig->SetParameter(1, omega_mass);

    // width
    sig->SetParameter(2, expected_width);



    TF1* bg = new TF1("bg", "pol2", r_min, r_max);
    bg->SetLineColor(kBlue);

    bg->SetParameter(0,0);
    bg->SetParName(0, "BG p_{0}");
    bg->SetParameter(1,0);
    bg->SetParName(1, "BG p_{1}");
    bg->SetParameter(2,0);
    bg->SetParName(2, "BG p_{2}");

    TFSum::FitRanged(h, bg, 650, 730, 830, 900);


    TFSum* sum = new TFSum("sum", sig, bg, r_min, r_max);
    sum->SetNpx(npx);


    //TCanvas* c = new TCanvas();
    //c->SetTitle(Form("Fit to %s", h->GetName()));
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

    cout << "Mass offset = " << peak_pos - omega_mass << " MeV\n";
    cout << "Sig/BG      = " << sig_to_bg << "\n";
    cout << "Sig         = " << sig_area << endl;

    // TODO: choose a position. Positions for TLatex are histogram coordinates.
    TLatex* label = new TLatex(r_min, h->GetMaximum(),Form("Signal content = %lf", sig_area));
    label->Draw();

    return FitResult(peak_pos, sig_area, sig->GetParameter(2), 0);

}

Omega::FitResult Omega::FitHistCrystalBall(TH1 *h, const bool fixOmegaMass, const double r_min, const double r_max) {

    const int    npx   = 500;
    const double omega_mass     = ParticleTypeDatabase::Omega.Mass();
    const double expected_width =  15.0;

    TF1* sig = ant::Math::CrystalBall();
    sig->SetLineColor(kGreen);
    sig->SetNpx(npx);

    constexpr auto iAlpha=0;
    constexpr auto iN=1;
    constexpr auto iSigma=1;
    constexpr auto iMass=3;

    // alpha
    sig->SetParameter(iAlpha, 1.0);

    // N
    sig->SetParameter(iN, 1.0);

    // position
    if(fixOmegaMass)
        sig->FixParameter(iMass, omega_mass);
    else
        sig->SetParameter(iMass, omega_mass);

    // width
    sig->SetParameter(2, expected_width);



    TF1* bg = new TF1("bg", "pol2", r_min, r_max);
    bg->SetLineColor(kBlue);

    bg->SetParameter(0,0);
    bg->SetParName(0, "BG p_{0}");
    bg->SetParameter(1,0);
    bg->SetParName(1, "BG p_{1}");
    bg->SetParameter(2,0);
    bg->SetParName(2, "BG p_{2}");

    TFSum::FitRanged(h, bg, 650, 730, 830, 900);


    TFSum* sum = new TFSum("sum", sig, bg, r_min, r_max);
    sum->SetNpx(npx);


    //TCanvas* c = new TCanvas();
    //c->SetTitle(Form("Fit to %s", h->GetName()));
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

    cout << "Mass offset = " << peak_pos - omega_mass << " MeV\n";
    cout << "Sig/BG      = " << sig_to_bg << "\n";
    cout << "Sig         = " << sig_area << endl;

    // TODO: choose a position. Positions for TLatex are histogram coordinates.
    TLatex* label = new TLatex(r_min, h->GetMaximum(),Form("Signal content = %lf", sig_area));
    label->Draw();

    return FitResult(peak_pos, sig_area, sig->GetParameter(2), 0);

}

TF1 *Omega::CB()
{
    return  ant::Math::CrystalBall();
}

struct FitDrawable : ant::root_drawable_traits {
    TH1* h;
    vector<Omega::FitResult>& res;

    FitDrawable(TH1* hist, vector<Omega::FitResult>& r): h(hist), res(r) {}

    void Draw(const string&) const override {
        res.emplace_back(Omega::FitHist(h));
    }
};

void Omega::FitOmegaPeak(const bool fixOmegaMass, const double r_min, const double r_max) {

    const char*  hist_name = "ggg_IM";


    TH1* h = NULL;
    gDirectory->GetObject(hist_name, h);

    if(!h) {
        cerr << "Can't find histogram" << endl;
        return;
    }

    FitHist(h, fixOmegaMass, r_min, r_max);

}

void Omega::FitOmega2D(TH2 *h) {

    canvas c;
    vector<Omega::FitResult> res;
    vector<double> bins;

    for(int i=2; i<h->GetNbinsY(); ++i) {
        TH1* p = h->ProjectionX(Form("%s_%d", h->GetName(), i),i,i+1);
        double bmin = h->GetYaxis()->GetBinLowEdge(i);
        double bmax = h->GetYaxis()->GetBinUpEdge(i);
        p->SetTitle(Form("%f - %f", bmin, bmax));
        c << FitDrawable(p,res);
        bins.push_back(h->GetYaxis()->GetBinCenter(i));
    }
    c << endc;

    TGraph* pos   = new TGraph(res.size());
    pos->SetTitle("#omega peak pos");

    TGraph* width = new TGraph(res.size());
    width->SetTitle("#omega peak width");

    TGraph* sig = new TGraph(res.size());
    sig->SetTitle("Signal");

    for(unsigned  i=0;i<res.size();++i) {
        pos->SetPoint(i,bins.at(i),res.at(i).pos);
        width->SetPoint(i,bins.at(i),res.at(i).width);
        sig->SetPoint(i,bins.at(i),res.at(i).N);
    }

    canvas() << pos << width << sig << endc;
}
