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
#include "analysis/plot/HistogramFactory.h"
#include "base/math_functions/CrystalBall.h"
#include "base/math_functions/AsymGaus.h"

#include "base/std_ext/memory.h"
#include "analysis/utils/MCWeighting.h"

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

struct SigFctWrapper {
    SigFctWrapper(TF1* f): sig(f) {
        sig->SetLineColor(kGreen);
        sig->SetNpx(500);
    }
    TF1* sig = nullptr;
    virtual double getPos() const =0;
    virtual double getWidth() const =0;
    virtual ~SigFctWrapper() = default;
};

struct SigFctGaus : SigFctWrapper {
    SigFctGaus(const double r_min, const double r_max, const double height, const double pos, const double width, const bool fixmass): SigFctWrapper(new TF1("sig", "gaus", r_min, r_max))
    {
        // height
        sig->SetParameter(0, 0.5 * height);

        // position
        if(fixmass)
            sig->FixParameter(1, pos);
        else
            sig->SetParameter(1, pos);

        // width
        sig->SetParameter(2, width);
    }

    double getPos() const override {
        return sig->GetParameter(1);
    }

    double getWidth() const override {
        return sig->GetParameter(2);
    }
};

struct SigFctCrystalBall : SigFctWrapper {

    static constexpr auto iAlpha=0;
    static constexpr auto iN=1;
    static constexpr auto iSigma=2;
    static constexpr auto iMass=3;
    static constexpr auto iHeight=4;

    SigFctCrystalBall(const double r_min, const double r_max, const double height, const double pos, const double width, const bool fixmass): SigFctWrapper(ant::Math::CrystalBall())
    {
        sig->SetRange(r_min,r_max);

        // alpha
        sig->SetParameter(iAlpha, 1.0);

        // N
        sig->SetParameter(iN, 1.0);

        // position
        if(fixmass)
            sig->FixParameter(iMass, pos);
        else
            sig->SetParameter(iMass, pos);

        // width
        sig->SetParameter(iSigma, width);

        sig->SetParameter(iHeight, height);
    }

    double getPos() const override {
        return sig->GetParameter(iMass);
    }

    double getWidth() const override {
        return sig->GetParameter(iSigma);
    }
};

Omega::FitResult Omega::FitHist(TH1 *h, const double omega_mass_, const bool fixOmegaMass, const FCT_t fct, const double r_min, const double r_max) {

    const int    npx   = 500;
    const double omega_mass     = omega_mass_ <= 0.0 ? ParticleTypeDatabase::Omega.Mass() : omega_mass_;
    const double expected_width =  15.0;

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

    const auto maxbin = h->GetMaximumBin();
    const auto height = h->GetBinContent(maxbin) - bg->Eval(h->GetBinCenter(maxbin));

    auto sig = [&]  () -> std::unique_ptr<SigFctWrapper> {
        switch (fct) {
        case eGAUS:
            return std_ext::make_unique<SigFctGaus>(r_min, r_max, height, omega_mass, expected_width, fixOmegaMass);
        case eCrystalBall:
            return std_ext::make_unique<SigFctCrystalBall>(r_min, r_max, height, omega_mass, expected_width, fixOmegaMass);
        }
        return nullptr;
    }();

    if(!sig) {
        cerr << "invalid singal function" << endl;
        return {};
    }





    TFSum* sum = new TFSum("sum", sig->sig, bg, r_min, r_max);
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
    const double peak_pos  = sig->getPos();
    const double chi2dnf = sum->Function()->GetChisquare()/sum->Function()->GetNDF();

    cout << "Mass offset = " << peak_pos - omega_mass << " MeV\n";
    cout << "Sig/BG      = " << sig_to_bg << "\n";
    cout << "Sig         = " << sig_area << endl;
    cout << "Chi2/dof    = " << chi2dnf << endl;

    // TODO: choose a position. Positions for TLatex are histogram coordinates.
    TLatex* label = new TLatex(r_min, h->GetMaximum(),Form("Signal content = %lf", sig_area));
    label->Draw();

    return FitResult(peak_pos, sig_area, sig->getWidth(), chi2dnf);

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

    FitHist(h, -1.0, fixOmegaMass, eGAUS, r_min, r_max);

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
#include "base/std_ext/memory.h"
#include "TRandom3.h"
TH2D* SampleDiffXsection(const ParticleTypeDatabase::Type& p)
{
    auto h = new TH2D("","Omega diff xsec bin sample prob", 470,1420,1580, 200,-1,1);
    analysis::HistogramFactory hf("a");
    unique_ptr<analysis::utils::MCWeighting> mcw = [&p,&hf] () {
        if(p==ParticleTypeDatabase::Omega)
            return std_ext::make_unique<analysis::utils::MCWeighting>(hf, analysis::utils::MCWeighting::Omega);
        if(p==ParticleTypeDatabase::Pi0)
            return std_ext::make_unique<analysis::utils::MCWeighting>(hf, analysis::utils::MCWeighting::Pi0);
        if(p==ParticleTypeDatabase::Eta)
            return std_ext::make_unique<analysis::utils::MCWeighting>(hf, analysis::utils::MCWeighting::Eta);
        return unique_ptr<analysis::utils::MCWeighting>(nullptr);
    }();

    if(!mcw) {
        cerr << "ERROR" << endl;
        return nullptr;
    }
//    for(Long64_t i=0;i<n;++i) {
//        const auto E = rnd.Uniform(1420,1580);
//        const auto cost = rnd.Uniform(-1,1);
//        const auto N = mcw->GetN(E,cost);
//        h->Fill(E,cost,N);
//    }

    for(int x=1;x<h->GetNbinsX(); ++x) {
        for(int y=1;y<h->GetNbinsY(); ++y) {
            h->SetBinContent(x,y,mcw->GetN(h->GetXaxis()->GetBinCenter(x), h->GetYaxis()->GetBinCenter(y)));
        }
    }

    return h;

}

TH2D* Omega::SampleDiffXsectionOmega()
{
    auto h = SampleDiffXsection(ParticleTypeDatabase::Omega);
    h->Scale(1.0/6.91090694237991715e-02);
    return h;
}

TH2D* Omega::SampleDiffXsectionEta()
{
    return SampleDiffXsection( ParticleTypeDatabase::Eta);
}
TH2D* Omega::SampleDiffXsectionPi0()
{
    auto h = SampleDiffXsection(ParticleTypeDatabase::Pi0);
    h->Scale(1.0/4.90141090503409238e+00);
    return h;
}

TH1D *Omega::getSignalYield(const string &meson)
{
    const auto n = 5;
    TH1D* yield = new TH1D("","",n,-1.0,1.0);
    for(int i=0;i<n; ++i) {
        TH1* h = nullptr;
        gDirectory->GetObject(Form("OmegaEtaG_Plot/Prob+mm/%sHyp/cosT_%d/h/Sum_MC/%shyp_omega",meson.c_str(), i, meson.c_str()), h);
        if(h) {
            yield->SetBinContent(i+1,h->Integral());
        }
    }
    return yield;
}
