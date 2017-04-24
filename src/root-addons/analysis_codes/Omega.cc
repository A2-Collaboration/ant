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
#include "analysis/plot/RootDraw.h"
#include "TGraph.h"
#include "base/ParticleType.h"

#include "base/math_functions/CrystalBall.h"
#include "base/math_functions/AsymGaus.h"

using namespace ant;
using namespace std;


inline double nCB_eval(double x, double A, double alpha, double n, double sigma, double mean, double hight) noexcept {
    return A * ant::math::CrystalBall::Eval(x,alpha,n,sigma,mean,hight);
}


inline double nCB_pol3_eval(double x, double A, double alpha, double n, double sigma, double mean, double hight, double p0, double p1, double p2, double p3) noexcept {
    auto v = nCB_eval(x,A,alpha,n,sigma,mean, hight) + p0;
    v += p1*x;
    x *= x;
    v += p2*x;
    x *= x;
    v += p3*x;
    return v;
}

inline double nCB_pol3_eval_ROOT(double* x, double* p) {
    return nCB_pol3_eval(x[0], p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9]);
}

inline double Pol_Eval_ROOT(double* x, double* p, int n) {
    double res = p[0];
    double xx = x[0];
    for(int i=1;i<n;++i){
        res += p[i] *xx;
        xx*=xx;
    }
    return res;
}

inline double AsymGaus_pol3_eval_ROOT(double* x, double* p) {
    return ant::math::AsymGaus::Eval_ROOT(x, p) + Pol_Eval_ROOT(x,(p+4),3);
}


//template<typename... Doubles>
//inline double eval(double x, Doubles... doubles) {
//    constexpr auto nDoubles = sizeof...(Doubles);

//}

//template<class Ch, class Tr, class Tuple, std::size_t... Is>
//void eval_impl(double* x, double* p, std::index_sequence<Is...>)
//{
//    (void)swallow{0, (void(os << (Is == 0 ? "" : ", ") << std::get<Is>(t)), 0)...};
//}


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

Omega::FitResult Omega::FitHist(TH1 *h, const double omega_mass_, const bool fixOmegaMass, const double r_min, const double r_max) {

    const int    npx   = 500;
    const double omega_mass     = omega_mass_ <= 0.0 ? ParticleTypeDatabase::Omega.Mass() : omega_mass_;
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

    const auto peak_range = interval<double>::CenterWidth(omega_mass, 4*expected_width);

    TFSum::FitRanged(h, bg, r_min, peak_range.Start(), peak_range.Stop(), r_max);


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
    constexpr auto iSigma=2;
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
    sig->SetParameter(iSigma, expected_width);



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

//TF1 *Omega::nCB()
//{
//    auto f = new TF1("", nCB_eval, -10,10, 5);
//    f->SetParNames("A","alpha","n","sigma","mean");
//    f->SetParameter(0,4.0);
//    f->SetParameter(1,1.0);
//    f->SetParameter(2,1.0);
//    f->SetParameter(3,1.0);
//    f->SetParameter(4,1.0);

//    return f;
//}

TF1 *Omega::nCB_pol3()
{
    auto f = new TF1("", nCB_pol3_eval_ROOT, -10,10, 5+4);
    f->SetParNames("A","alpha","n","sigma","mean","p0","p1","p2","p3");
    f->SetParameter(0,3.29023e+04);
    f->SetParameter(1,-1.54640e+00);
    f->SetParameter(2,3.00000e+00);
    f->SetParameter(3,1.67039e+01);
    f->SetParameter(4,782);
    f->SetParameter(5,0);
    f->SetParameter(6,0);
    f->SetParameter(7,0);
    f->SetParameter(8,0);

    return f;
}

TF1 *Omega::AsymGaus_pol3()
{
    auto f = new TF1("", AsymGaus_pol3_eval_ROOT, -10,10, 4+4);
    f->SetParNames("A","mean","s1","s2","p0","p1","p2","p3");
    f->SetParameter(0,3.29023e+04);
    f->SetParameter(1,782);
    f->SetParameter(2,1);
    f->SetParameter(3,1);
    f->SetParameter(4,0);
    f->SetParameter(5,0);
    f->SetParameter(6,0);
    f->SetParameter(7,0);
    return f;
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
