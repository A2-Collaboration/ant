#pragma once

#include <vector>
#include <algorithm>
#include <iostream>

#include "APLCON.hpp"


#include "base/interval.h"
#include "base/Logger.h"
#include "base/ParticleType.h"
#include "base/WrapTFile.h"


#include "base/std_ext/math.h"

#include "analysis/utils/ValError.h"

#include "TGraphErrors.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"

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



class TGraphErrors;
class TGraph;
class TF1;

using namespace ant;
using namespace std;
using namespace RooFit;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::utils;

namespace SIGMA {

struct BR{
    static constexpr auto is_K0S     = 0.5;
    static constexpr auto Sigma_Pi0p = 0.5157;
    static constexpr auto K0S_Pi0Pi0 = 0.3069;
    static constexpr auto Pi0_gg     = 0.98823;
};

struct functions{
    static ValError MakeDataPoint(const RooRealVar& value)
    {
        return ValError(value.getValV(),value.getError());
    }

    template <typename T>
    static T* getHist(const WrapTFileInput& f, const string& hpath) {
        T* h = nullptr;
        if(!f.GetObject(hpath, h)) {
            LOG(FATAL) << "Cannot find " << hpath;
        };
        return h;
    }

    static vector<TH2D*> getHists(WrapTFileInput& f, const string& hpath,const string& sndpath, const string& histname, const unsigned nbins)
    {
        vector<TH2D*> hs;
        const auto paths = [&hpath,&sndpath,&histname,nbins]()
        {
            vector<string> ret;
            for (auto egbin = 0u ; egbin < nbins; ++egbin)
            {
                const string histp = std_ext::formatter() << hpath << egbin << "/" << sndpath << histname;
                ret.push_back(histp);
            }
            return ret;
        }();

        transform(paths.begin(),paths.end(),
                  back_inserter(hs),[&f](const string& p){return getHist<TH2D>(f,p);});

        return hs;
    }

    static vector<ValError> getSeenMCsums(const WrapTFileInput& f,
                                          const string& hpath, const string& prefix,
                                          const unsigned nEgbins)
    {
        vector<ValError> sums;
        sums.reserve(nEgbins);
        for ( auto egbin = 0u; egbin < nEgbins; ++egbin)
        {
            const string histp = std_ext::formatter() << hpath << prefix << egbin;
            auto h = getHist<TH1D>(f,histp);
            sums.push_back(ValError::Statistical(h->Integral()));
        }
        return  sums;
    }
};

struct tools {
    static size_t FillGraphErrors(TGraphErrors* graph, const double x, const double y, const double xerr, const double yerr)
    {
        auto N = graph->GetN();
        graph->SetPoint(N,x,y);
        graph->SetPointError(N,xerr,yerr);
        return graph->GetN();
    }

    static ValError fitHist (TH2D* histDalitz)
    {
        RooRealVar mppi0("mppi0","m(p #pi^{0})",1150, 1230) ;
        RooRealVar mpi0pi0("mpi0pi0","m(#pi^{0} #pi^{0})", 420, 560) ;

        RooDataHist dh ("dh", "dh", RooArgSet(mppi0, mpi0pi0), histDalitz);


        RooRealVar mean_Sigma("mean_Sigma","m_{#Sigma}",1188, 1170, 1220);
        RooRealVar sigma_Sigma("sigma_Sigma","#sigma_{#Sigma}",10, 0, 50);
        RooRealVar alpha_Sigma("alpha_Sigma","#alpha_{#Sigma}",2);
        RooRealVar n_Sigma("n_Sigma","n_Sigma",1);
        RooCBShape *theSigma = new RooCBShape("theSigma","crystal ball PDF",
                                              mppi0,mean_Sigma,sigma_Sigma,alpha_Sigma,n_Sigma);

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

        RooFitResult *fitRes = fullPdf->fitTo (dh, SumW2Error(kTRUE),
                                               Minimizer("Minuit"), PrintLevel(1),
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



        // Create and fill ROOT 2D histogram (8x8x8 bins) with contents of dataset
        TH1* hh_data3 = dh.createHistogram("hh_data3",
                                           mppi0,Binning(24),
                                           YVar(mpi0pi0, Binning(21))) ;


        // Create and fill ROOT 2D histogram (20x20x20 bins) with sampling of pdf
        TH1* hh_pdf3 = fullPdf->createHistogram("hh_model3",
                                                mppi0,Binning(100),
                                                YVar(mpi0pi0,Binning(100))) ;
        hh_pdf3->SetFillColor(kBlue) ;

        if (fitRes->status() == 0)
            return {nSig.getVal(),nSig.getError()};
        return {NaN,NaN};
    }

//    static ValError GetLumi(TH1D* hlumi, const IntervalD& egRange)
//    {

//        return ValError(1,1);
//    }

};

struct CrossSectionDataPoint
{
    double Egamma;
    ValError Data;
    string Unit;

    CrossSectionDataPoint(const double eGamma, const ValError data,
                          const string& unit = "#mub"):
        Egamma(eGamma), Data(data), Unit(unit) {}
    CrossSectionDataPoint(const double eGamma, const RooRealVar& value,
                          const string& unit = "#mub"):
        Egamma(eGamma), Data(functions::MakeDataPoint(value)), Unit(unit) {}

    double W() const {return sqrt(
                    sqr(Egamma + ParticleTypeDatabase::Proton.Mass()) - sqr(Egamma));}

};

struct LumiFitter_t
{

    struct Value_t {

        constexpr Value_t(double v, double s) : V_S_P{v, s, 0}{}

        std::tuple<double,double,double> V_S_P; // short for Value, Sigma, Pull

        bool Fixed = false;
        bool Poisson = false;

        // fitter is able to access fields to fit them (last one pull is just written, actually)
        template<size_t N>
        std::tuple<double&> linkFitter() noexcept {
            return std::tie(std::get<N>(V_S_P));
        }

        // user can provide settings for each variable
        template<size_t innerIdx>
        APLCON::Variable_Settings_t getFitterSettings(size_t outerIdx) const noexcept {
            (void)outerIdx; // unused, provided to user method for completeness
            APLCON::Variable_Settings_t settings;
            if(Fixed)
                settings.StepSize = 0;
            if(Poisson)
                settings.Distribution = APLCON::Distribution_t::Poissonian;
            return settings;
        }

        // the following methods are just user convenience (not required by fitter)

        friend std::ostream& operator<<(std::ostream& s, const Value_t& o) {
            return s << "(" << o.Value() << "," << o.Sigma() << "," << o.Pull() << ")";
        }

        double Value() const {
            return std::get<APLCON::ValueIdx>(V_S_P);
        }

        double Sigma() const {
            return std::get<APLCON::SigmaIdx>(V_S_P);
        }

        double Pull() const {
            return std::get<APLCON::PullIdx>(V_S_P);
        }

        operator double() const {
            return std::get<APLCON::ValueIdx>(V_S_P);
        }
    };

    APLCON::Fitter<Value_t,Value_t,Value_t,Value_t,Value_t,Value_t,
                   std::vector<Value_t>, std::vector<Value_t>> fitter;

    LumiFitter_t() = default;

    struct result_t
    {
        Value_t Intercept{0,0};
        Value_t c1{0,0};
        Value_t c2{0,0};
        Value_t c3{0,0};
        Value_t c4{0,0};
        Value_t c5{0,0};

        std::vector<Value_t> Fitted_Egs;
        std::vector<Value_t> Fitted_Lumis;
        APLCON::Result_t FitStats;
    };
    LumiFitter_t::result_t DoFit(TH1D* histLumi);

    static TF1* makeROOTfunction(const LumiFitter_t::result_t& result, const double& xmin, const double& xmax)
    {
        auto fu = [result](double* xs, double*)
        {
            return result.Intercept
                    + result.c1 * xs[0]
                    + result.c2 * pow(xs[0],2)
                    + result.c3 * pow(xs[0],3)
                    + result.c4 * pow(xs[0],4)
                    + result.c5 * pow(xs[0],5);
        };
        return new TF1("fit",fu,xmin,xmax,0);
    }
};



}
