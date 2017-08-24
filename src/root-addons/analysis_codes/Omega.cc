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
#include "base/std_ext/math.h"

#include "base/std_ext/memory.h"
#include "analysis/utils/MCWeighting.h"
#include "root-addons/analysis_codes/hstack.h"
#include "TNtupleD.h"
#include "TGraphErrors.h"
#include "TGaxis.h"
#include "base/TH_ext.h"
#include "base/PhysicsMath.h"

using namespace ant;
using namespace ant::std_ext;
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

    // positionCanvas_1_n2
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


void Omega::SaveStacks(const string &path_spec, const string& fname, const int start, const int stop)
{
    TCanvas* c = new TCanvas();
    c->SetCanvasSize(600,600);
    c->SetTicks(1,1);

    for(int i=start; i<=stop; ++i) {
        const string path = Form(path_spec.c_str(), i);
        hstack* s = nullptr;
        gDirectory->GetObject(path.c_str(), s);
        if(!s) {
            cerr << path << " not found" << endl;
        } else {
            c->cd();
            auto l = new TLatex( 650.0, s->GetMaximum()*0.9, Form("cos(#theta) = %f", -0.9+0.2*i));
            s->Draw();
            l->Draw();
            c->SaveAs(Form(fname.c_str(),i));
        }
    }
}

inline TGraph* makeGraph(const vector<vec2>& data) {

    auto g = new TGraph(int(data.size()));
    for(int i=0; i<int(data.size()); ++i) {
        const auto& p = data.at(size_t(i));
        g->SetPoint(i, p.x, p.y);
    }
    return g;
}

TGraph* getRefGraph(const double& W, const vector<pair<interval<double>,vector<vec2>>>& data) {
    for(const auto& rg : data) {
        if(rg.first.Contains(W)) {
            auto g = makeGraph(rg.second);
            g->SetTitle(Form("W=%.1f", rg.first.Center()));
            return g;
        }
    }

    return nullptr;
}



void Omega::PlotFitted(const string &file)
{
    auto t = new TNtupleD("omegafitdata","","cosT:Emin:Emax:Ecenter:Nsig:dNsig:Nbkg:dNbkg:RecEff:dRecEff:Ncrr:dNcorr:nMC:nMCInput:sigma:dsigma:chi2dof:mshift:gauss:c:p:f");
    //auto t = new TNtupleD("omegafitdata","","cosT:Emin:Emax:Ecenter:Nsig:dNsig:Nbkg:dNbkg:RecEff:dRecEff:Ncrr:dNcorr:DataToMC:sigma:dsigma");
    t->ReadFile(file.c_str());

    double cosT, Ecenter,Emin,Emax,RecEff,Ncorr,nMCInput,sigma,dSigma,mshift,gauss;
    t->SetBranchAddress("cosT",&cosT);
    t->SetBranchAddress("Emin",&Emin);
    t->SetBranchAddress("Emax",&Emax);
    t->SetBranchAddress("Ecenter",&Ecenter);
    t->SetBranchAddress("RecEff",&RecEff);
    t->SetBranchAddress("Ncrr",&Ncorr);
    t->SetBranchAddress("sigma",&sigma);
    t->SetBranchAddress("dsigma",&dSigma);
    t->SetBranchAddress("nMCInput",&nMCInput);
    t->SetBranchAddress("mshift",&mshift);
    t->SetBranchAddress("gauss",&gauss);

    t->GetEntry(0);
    cout << Ecenter << endl;

    TH2D* s2d        = new TH2D("sigma2d" ,"cross section 2d",10,-1,1,12,1420,1580);
    TH2D* h_nMCInput = new TH2D("nMCInput","n MC Input 2d"   ,10,-1,1,12,1420,1580);
    TH2D* h_Ncorr    = new TH2D("Ncorr",   "Corrected Number of events"   ,10,-1,1,12,1420,1580);
    TH2D* h_RecEff   = new TH2D("RecEff",  "Reconstruction Efficiency "   ,10,-1,1,12,1420,1580);
    TH2D* h_mshift   = new TH2D("mshift",  "Omega Mass Offset"            ,10,-1,1,12,1420,1580);
    TH2D* h_gauss    = new TH2D("gauss",  "Gauss Smearing"                ,10,-1,1,12,1420,1580);
    const auto setBinXY = [] (TH2* h, const double x, const double y, const double v) {
      const auto bin = h->FindBin(x,y);
      h->SetBinContent(bin,v);
    };

    map<double,TH1D*> graphs;

    for(int i=0;i<t->GetEntries(); ++i) {
        t->GetEntry(i);

        auto g = [&graphs] (const double Ec) {
            auto ge = graphs.find(Ec);
            if(ge!=graphs.end()) {
                return ge->second;
            } else {
                auto ng = new TH1D("","", 10, -1, 1);
                ng->GetYaxis()->SetRangeUser(0,10);
                ng->SetStats(false);
                ng->SetTitle("This work");
                ng->SetBit(TH1::kNoTitle);
                ng->SetXTitle("cos(#theta)_{cm}^{#omega}");
                ng->SetYTitle("#frac{d#sigma}{d cos(#theta)_{cm}} [#mub]");
                ng->SetLineColor(kBlack);
                ng->GetListOfFunctions()->Add(new TLatex(-.9,9,  Form("W=%.1f MeV", math::W(Ec,ParticleTypeDatabase::Proton))));
                ng->GetListOfFunctions()->Add(new TLatex(-.9,8.3,Form("E_{#gamma}=%.1f MeV", Ec)));

                ng->GetXaxis()->SetLabelSize(0.05f);
                ng->GetXaxis()->SetTitleSize(0.05f);
                ng->GetXaxis()->SetTitleOffset(0.95f);
                ng->GetXaxis()->SetNdivisions(509);

                ng->GetYaxis()->SetLabelSize(0.05f);
                ng->GetYaxis()->SetTitleSize(0.05f);
                ng->GetYaxis()->SetTitleOffset(0.95f);
                ng->GetYaxis()->SetNdivisions(509);
                graphs[Ec] = ng;
                return ng;
            }
        }(Ecenter);

        const auto bin = g->FindBin(cosT);
        g->SetBinContent(bin,sigma);
        g->SetBinError(bin,dSigma);

        setBinXY(s2d,cosT,Ecenter,sigma);
        setBinXY(h_nMCInput,cosT,Ecenter, nMCInput);
        setBinXY(h_Ncorr,cosT,Ecenter, Ncorr);
        setBinXY(h_RecEff,cosT,Ecenter,RecEff);
        setBinXY(h_mshift,cosT,Ecenter,mshift);
        setBinXY(h_gauss,cosT,Ecenter,gauss);

        cout << Ecenter << " " << cosT << " " << nMCInput << sigma << endl;
    }

    auto save_c = new TCanvas();
    save_c->SetCanvasSize(400,400);
    save_c->SetTicks(1,1);
    save_c->SetMargin(.15f,.02f,.15f,.02f);

    const vector<pair<interval<double>,vector<vec2>>> clasData = {
        {interval<double>::CenterWidth(1945,10),{
             { -0.7545123988600251 , 1.9731385331018423  },
             { -0.6454286369443434 , 1.9770694794771835   },
             { -0.5512169554820314 , 2.160644675205557   },
             { -0.45674321092803094 , 2.2000851705047992   },
             { -0.34071477708258247 , 2.384446555508241   },
             { -0.25351328332295964 , 2.423624987715794   },
             { -0.15163625642873502 , 2.391260195892162   },
             { -0.04982474530743275 , 2.3949290791758138   },
             { 0.052314344678481284 , 2.218429586923051   },
             { 0.15412585579978444 , 2.222098470206703   },
             { 0.2562649457856985 , 2.0455989779539436   },
             { 0.3578799095882337 , 2.1573688865594356   },
             { 0.4597569364824592 , 2.1250040947358073   },
             { 0.5531824286697029 , 2.7409833917515662   },
             { 0.653487077013791 , 3.573426802502702   },
             { 0.7454712221967448 , 4.982146951878665   },
         }},
        {interval<double>::CenterWidth(1935,10),{
             { -0.7508688376882482 , 2.0311494400823804   },
             { -0.641845797399923 , 1.9938215986613486   },
             { -0.5396447419230279 , 2.2044021109537972   },
             { -0.43808726991890845 , 2.0610117132192087   },
             { -0.33575749774745756 , 2.3423864075170613   },
             { -0.24089329386021374 , 2.517698545501357   },
             { -0.13933582185609517 , 2.3743081477667722   },
             { -0.03758527481014262 , 2.337109023040295   },
             { 0.06410091388853134 , 2.2645128073111174   },
             { 0.1583858926502768 , 2.12125112627108   },
             { 0.2600720813489512 , 2.0486549105419023   },
             { 0.3616295533530698 , 1.9052645128073173   },
             { 0.463894967177243 , 2.1512421161024626   },
             { 0.5665465310850819 , 2.6096022654138302   },
             { 0.6696486034238642 , 3.3157420517441203   },
             { 0.7661217659930495 , 4.375981464795991   },
         }},
        {interval<double>::CenterWidth(1925,10),{
             { -0.7497205234431513 , 2.0720720720720713 },
             { -0.6473334648517132 , 2.18018018018018 },
             { -0.5450779246399686 , 2.2162162162162176 },
             { -0.4428223844282241 , 2.2522522522522515 },
             { -0.3399092523180116 , 2.648648648648649 },
             { -0.2452817781284936 , 2.5045045045045065 },
             { -0.14289471953705546 , 2.6126126126126117 },
             { -0.048398763727231175 , 2.396396396396394 },
             { 0.06115604655750606 , 2.4324324324324316 },
             { 0.15578352074702417 , 2.288288288288289 },
             { 0.2577760241993816 , 2.18018018018018 },
             { 0.35246925757874625 , 2.0720720720720713 },
             { 0.4548563161701846 , 2.18018018018018 },
             { 0.5653975143026233 , 2.7567567567567544 },
             { 0.6619320049976982 , 3.6576576576576585 },
             { 0.7577431446044582 , 4.162162162162161 },
         }},
        {interval<double>::CenterWidth(1915,10),{
             { -0.7410071942446042 , 2.2026279802604165 },
             { -0.6474820143884887 , 2.3033474047208493 },
             { -0.5467625899280577 , 2.4372435935549035 },
             { -0.4460431654676258 , 2.5050240799096244 },
             { -0.34532374100719476 , 2.7050359712230154 },
             { -0.25179856115107935 , 2.77269754444378 },
             { -0.15107913669064743 , 2.509899518401806 },
             { -0.057553956834532904 , 2.544503240382898 },
             { 0.050359712230215514 , 2.4140555324335544 },
             { 0.14388489208633093 , 2.448659254414647 },
             { 0.25179856115107935 , 2.2851536952256346 },
             { 0.3525179856115104 , 2.2207027766216747 },
             { 0.4532374100719432 , 2.4537725191747377 },
             { 0.5467625899280586 , 2.9511861585112023 },
             { 0.6474820143884887 , 3.8123550746179866 },
             { 0.7482014388489215 , 4.442119032047083 },
         }},
        {interval<double>::CenterWidth(1905,10),{
             { -0.7515603592632063 , 2.1975338712132775 },
             { -0.6499314964225915 , 2.2966661592327604 },
             { -0.5624904856142496 , 2.564104125437666 },
             { -0.44679555487897726 , 2.42770589130766 },
             { -0.3520474958136708 , 2.728634495356985 },
             { -0.2435378291977477 , 2.6259704673466295 },
             { -0.15670573907748508 , 2.5572842137311653 },
             { -0.04825696453037054 , 2.4210077637387784 },
             { 0.05349368244786046 , 2.587364895722338 },
             { 0.14738925254985435 , 2.4177195920231433 },
             { 0.24840919470238898 , 2.18072766022226 },
             { 0.3497944892677731 , 2.1454102603135965 },
             { 0.45178870452123565 , 2.4462170802253027 },
             { 0.5542091642563554 , 2.982310854011267 },
             { 0.6565687319226665 , 3.4847922058151966 },
             { 0.7597807885522903 , 4.4578474653676405 },

         }},
        {interval<double>::CenterWidth(1895,10),{
             { -0.6463529197854476 , 2.59395739264826 },
             { -0.5521076396254432 , 2.759053304645594 },
             { -0.45743810418497555 , 2.690808812388255 },
             { -0.3559804842570986 , 2.889117852056124 },
             { -0.24694687717809627 , 2.9206339586048067 },
             { -0.1451862177641745 , 2.952271280947908 },
             { -0.050153034940452645 , 2.6840207279008474 },
             { 0.051425800781842224 , 2.815661080638808 },
             { 0.15379253916785363 , 2.513954968332378 },
             { 0.24840146671111274 , 2.4790448195399897 },
             { 0.350283341919452 , 2.4440134549531827 },
             { 0.4520440013333742 , 2.475650777296284 },
             { 0.5456225946240796 , 3.007424467408107 },
             { 0.6533834358616932 , 3.7389617867208127 },
             { 0.746234734385891 , 4.670747598412074 },
         }},
        {interval<double>::CenterWidth(1885,10),{
             { -0.7483301012712779 , 2.504201680672267 },
             { -0.6531535691199558 , 2.470588235294116 },
             { -0.5500353987748946 , 2.7731092436974762 },
             { -0.4477175485578849 , 2.6386554621848717 },
             { -0.3447225043863701 , 2.873949579831933 },
             { -0.24985378766891386 , 2.672268907563023 },
             { -0.14722812201803814 , 2.705882352941174 },
             { -0.051559085172530494 , 2.9411764705882355 },
             { 0.05045094961061336 , 2.6386554621848717 },
             { 0.15283036291439633 , 2.5378151260504183 },
             { 0.2553329023917257 , 2.504201680672267 },
             { 0.35777387878228195 , 2.436974789915965 },
             { 0.4533197894542431 , 2.6050420168067205 },
             { 0.5566842121463969 , 3.0420168067226854 },
             { 0.6610336442269213 , 4.016806722689074 },
             { 0.7646443192661676 , 4.588235294117647 },
         }},
    };
    const vector<pair<interval<double>,vector<vec2>>> saphirData = {
        {interval<double>::CenterWidth(1885,10),{
             { -0.921497992419979 , 1.29652 },
             { -0.750909977860333 , 2.29384 },
             { -0.552253367856205 , 2.49331 },
             { -0.382003077038538 , 2.59304 },
             { -0.197643438778191 , 2.49331 },
             { -0.013021126496304 , 3.0917 },
             { 0.178355660625164 , 2.79251 },
             { 0.362602724304851 , 2.39358 },
             { 0.490262298772937 , 2.39358 },
             { 0.58981575293632 , 3.0917 },
             { 0.682464632819242 , 4.28849 },
             { 0.775188562422606 , 5.68475 },
             { 0.861420691207925 , 8.67672 },
             { 0.940448046831026 , 11.3695 },
         }},
        {interval<double>::CenterWidth(1898,10),{
             { -0.539007092198582 , 1.67211 },
             { -0.361702127659575 , 1.36407 },
             { -0.184397163120567 , 2.65175 },
             { -0.007092198581561 , 2.14425 },
             { 0.163120567375887 , 2.03603 },
             { 0.340425531914894 , 2.62558 },
             { 0.5177304964539 , 3.01567 },
             { 0.695035460992908 , 5.20094 },
             { 0.865248226950355 , 7.78549 },
         }},
        {interval<double>::CenterWidth(1910,10),{
             { -0.935943060498221 , 1.09128 },
             { -0.750889679715303 , 1.78573 },
             { -0.580071174377225 , 2.08336 },
             { -0.416370106761566 , 2.38098 },
             { -0.252669039145908 , 2.38098 },
             { -0.08896797153025 , 2.57939 },
             { 0.081850533807828 , 2.28177 },
             { 0.245551601423486 , 2.38098 },
             { 0.416370106761566 , 2.08336 },
             { 0.544483985765124 , 2.48019 },
             { 0.637010676156584 , 2.77781 },
             { 0.715302491103203 , 3.76988 },
             { 0.793594306049822 , 6.05166 },
             { 0.886120996441281 , 9.3255 },
             { 0.94306049822064 , 12.0041 },
         }},
        {interval<double>::CenterWidth(1934,10),{
             { -0.92164809186086 , 0.897591 },
             { -0.765357049044992 , 1.59572 },
             { -0.602086382228228 , 1.99465 },
             { -0.446095538294122 , 1.89492 },
             { -0.297046793500694 , 2.19411 },
             { -0.155052722428609 , 2.59304 },
             { 0.007917745506397 , 2.19411 },
             { 0.163871064580283 , 1.99465 },
             { 0.305715036211488 , 1.99465 },
             { 0.468985703028256 , 2.39358 },
             { 0.575293632031221 , 2.19411 },
             { 0.660962887913241 , 3.6901 },
             { 0.739690044654583 , 5.58501 },
             { 0.811249953093924 , 7.28046 },
             { 0.883823032759203 , 11.6687 },
             { 0.948515891778303 , 13.9625 },
         }},
        {interval<double>::CenterWidth(1959,10),{
             { -0.921708185053381 , 1.06381 },
             { -0.779359430604982 , 1.32977 },
             { -0.629893238434164 , 1.19679 },
             { -0.487544483985765 , 1.72869 },
             { -0.345195729537366 , 1.59572 },
             { -0.209964412811388 , 1.99465 },
             { -0.06761565836299 , 2.2606 },
             { 0.081850533807829 , 2.12762 },
             { 0.224199288256228 , 1.72869 },
             { 0.373665480427046 , 1.46274 },
             { 0.516014234875445 , 2.12762 },
             { 0.622775800711744 , 2.79251 },
             { 0.686832740213524 , 3.45739 },
             { 0.765124555160142 , 5.71799 },
             { 0.829181494661922 , 8.37752 },
             { 0.893238434163701 , 11.9679 },
             { 0.957295373665481 , 16.2231 },
         }},
    };

    canvas s("Sigma");

    auto save_multi_images = [] (const TCanvas* c, const char* basename) {
        c->SaveAs(Form("%s.pdf", basename));
        c->SaveAs(Form("%s.png", basename));
        c->SaveAs(Form("%s.root", basename));
    };

    int i=0;
    for(auto& e : graphs) {
        const auto w = math::W(e.first, ParticleTypeDatabase::Proton);

        auto clasGraph = getRefGraph(w, clasData);
        if(clasGraph) {
            clasGraph->SetTitle(Form("CLAS 2009, %s", clasGraph->GetTitle()));
            clasGraph->SetMarkerColor(kGreen);
            clasGraph->SetMarkerStyle(kOpenCircle);
            clasGraph->SetFillColor(kWhite);
            clasGraph->SetLineColor(kWhite);
        }
        auto saphirGraph = getRefGraph(w, saphirData);
        if(saphirGraph) {
            saphirGraph->SetTitle(Form("SAPHIR 2003, %s", saphirGraph->GetTitle()));
            saphirGraph->SetMarkerColor(kRed);
            saphirGraph->SetMarkerStyle(kOpenSquare);
            saphirGraph->SetFillColor(kWhite);
            saphirGraph->SetLineColor(kWhite);
        }

        s << e.second;
        save_c->cd();
        e.second->Draw();

        const auto drawif = [&] (TGraph* g) {
            if(g) {
                g->Draw("same P");
                s << samepad << drawoption("P") << g;
            }
        };

        drawif(clasGraph);
        drawif(saphirGraph);

        save_c->BuildLegend(0.175,0.5,0.6,0.75);
        save_multi_images(save_c, Form("cross_Sec_%d", i++));
    }
    s << endc;

    canvas() << drawoption("colz") << s2d << h_nMCInput << h_Ncorr << h_mshift << h_gauss << endc;

    delete t;
}

TH2D *Omega::ConvertMesonCountHistToCrossSec(const TH2D *)
{
    auto res = new TH2D("","",10,-1,1,12,1420,1580);
    return res;
}

TCanvas *Omega::Printhstack(const string &name)
{
    ant::hstack* s = nullptr;
    gDirectory->GetObject(name.c_str(), s);
    if(s)
        return Printhstack(s);
    return nullptr;
}

TCanvas *Omega::Printhstack(hstack *hs)
{
    auto c = new TCanvas();
    c->SetCanvasSize(600,600);
    c->SetTicks(1,1);
    c->SetMargin(0.15f,0.05f,0.15f,0.05f);

    TGaxis::SetMaxDigits(3);

    const auto setAxis = [] (TAxis* a) {
      a->SetLabelSize(0.05f);
      a->SetTitleSize(0.05f);
      a->SetTitleOffset(1.05f);
      a->SetNdivisions(505);
    };

    hs->Draw();

    setAxis(hs->GetXaxis());
    setAxis(hs->GetYaxis());

    return c;
}

void Omega::RoundBins(TH2 *h, int d)
{
    const double f = pow(10,d);
    for(int x=1;x<=h->GetNbinsX();++x)
        for(int y=1;y<=h->GetNbinsY();++y) {
            h->SetBinContent(x,y,round( h->GetBinContent(x,y) * f) / f );
        }
}

TH1D *Omega::HistOfBins(const TH2 *h2)
{
    return TH_ext::HistOfBins(h2);
}

std::string Omega::TH1ToLaTeX(const TH1 *h, const int p)
{
    return TH_ext::TH1ToLaTeX(h,p);
}
