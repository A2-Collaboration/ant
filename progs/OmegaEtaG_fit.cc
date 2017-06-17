#include "base/Logger.h"

#include "tclap/CmdLine.h"
#include "base/interval.h"
#include "base/WrapTFile.h"
#include "base/std_ext/string.h"
#include "base/std_ext/system.h"
#include "base/std_ext/memory.h"
#include "base/ParticleType.h"
#include "base/TH_ext.h"
#include "base/std_ext/string.h"
#include "base/math_functions/Linear.h"
#include "base/std_ext/math.h"
#include "analysis/plot/RootDraw.h"
#include "root-addons/analysis_codes/Math.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPaveText.h"

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
#include "RooPolynomial.h"
#include "TGraphErrors.h"


using namespace ant;
using namespace ant::std_ext;
using namespace std;
using namespace RooFit;

struct ValError {
    double v;
    double e;

    ValError(const double& value=NaN, const double& error=NaN):v(value),e(error) {}

    ValError(const RooRealVar& value): v(value.getValV()), e(value.getError()) {}

    ValError(const ValError&) = default;
    ValError(ValError&&) = default;
    ValError& operator=(const ValError&) = default;
    ValError& operator=(ValError&&) = default;

    friend ostream& operator<<(ostream& s, const ValError& v);

    ValError operator/ (const double d) const {
        return {v/d, e/d};
    }

    ValError operator/ (const ValError& o) const {
        return {
            v / o.v,
            sqrt(sqr(e / o.v) + sqr(v/sqr(o.v)*o.e))
        };
    }

    static ValError Statistical(const double v) {
        return {v, sqrt(v)};
    }
};

ostream& operator<<(ostream& s, const ValError& v) {
    s << v.v << " +/- " << v.e;
    return s;
}

struct FitOmegaPeak {
    static constexpr auto BRomegaPi0ggg = 0.081144;

    ValError vnsig = std_ext::NaN;
    ValError vnbkg = std_ext::NaN;
    int numParams = 0;
    double chi2ndf = std_ext::NaN;

    ValError rec_eff = std_ext::NaN;
    ValError vn_corr = std_ext::NaN;
    ValError sigmaOmega = std_ext::NaN;

    double DataToMc = std_ext::NaN;

    FitOmegaPeak() {}
    FitOmegaPeak(const TH1* data, const TH1* mc_shape, const double n_mc_input, const ValError& total_lumi, const double binwidth=1.0);
    FitOmegaPeak(const FitOmegaPeak&) = default;
    FitOmegaPeak(FitOmegaPeak&&) = default;
    FitOmegaPeak& operator=(const FitOmegaPeak&) = default;
    FitOmegaPeak& operator=(FitOmegaPeak&&) = default;

    friend ostream& operator<<(ostream& s, const FitOmegaPeak& f);
};

TH2D* extrudeX(const TH1* slice, const BinSettings& xbins, const string& newname="") {
    auto h = new TH2D(newname.c_str(),"", int(xbins.Bins()), xbins.Start(), xbins.Stop(), slice->GetNbinsX(), slice->GetXaxis()->GetXmin(), slice->GetXaxis()->GetXmax());
    const auto SetSlice = [] (TH2D* h, const int bin, const TH1* slicedata) {
        for(int i=0;i<=slicedata->GetNbinsX();++i) {
            h->SetBinContent(bin,i,slicedata->GetBinContent(i));
            h->SetBinError(bin,i,slicedata->GetBinError(i));
        }
    };

    for(int i=0;i<h->GetNbinsX();++i) {
        SetSlice(h,i,slice);
    }
    return h;
}

template <typename T>
T* getHist(WrapTFileInput& f, const string& hpath) {
    T* h;
    if(!f.GetObject(hpath, h)) {
        LOG(FATAL) << "Cannot find " << hpath;
    };
    return h;
};

using cosTbins_t = vector<tuple<double,double,FitOmegaPeak>>;

int main(int argc, char** argv) {
    SetupLogger();

    TCLAP::CmdLine cmd("OmegaEtaG_fit", ' ', "0.1");
    auto cmd_verbose   = cmd.add<TCLAP::ValueArg<int>>   ("v","verbose","Verbosity level (0..9)", false, 0,"int");
    auto cmd_data      = cmd.add<TCLAP::ValueArg<string>>("", "data","Data input",true,"","rootfile");
    auto cmd_lumi      = cmd.add<TCLAP::ValueArg<string>>("", "lumi","Lumi data",true,"","rootfile");
    auto cmd_mc        = cmd.add<TCLAP::ValueArg<string>>("", "mc","MC signal/reference input",true,"","rootfile");
    auto cmd_histpath  = cmd.add<TCLAP::ValueArg<string>>("", "histpath","Path for hists",false,"OmegaEtaG_Plot/n==4+Prob/pi0Hyp","path");
    auto cmd_histname  = cmd.add<TCLAP::ValueArg<string>>("", "histname","Name of hist",false,"ggg_IM","name");
    auto cmd_batchmode = cmd.add<TCLAP::MultiSwitchArg>  ("b","batch","Run in batch mode (no ROOT shell afterwards)",false);
    auto cmd_output    = cmd.add<TCLAP::ValueArg<string>>("o","output","Output file",false,"","filename");
    auto cmd_mc_inputn = cmd.add<TCLAP::ValueArg<string>>("", "mcinput","MC input file (generated nums)",true,"","filename");
    auto cmd_ITtest    = cmd.add<TCLAP::SwitchArg>  ("","iotest","Run input/output test",false);
    auto cmd_mode      = cmd.add<TCLAP::ValueArg<string>>("m","mode","Fit Mode: global, cosT, cosTE",true,"","mode");

    cmd.parse(argc, argv);
    if(cmd_verbose->isSet()) {
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());
    }

    const string datahist = cmd_ITtest->isSet() ? "/h/Sum_MC/" : "/h/Data/";
    const string refhist  = "/h/Ref/";


    WrapTFileInput input_lumi(cmd_lumi->getValue());

    WrapTFileInput input_mcnumbers(cmd_mc_inputn->getValue());

    auto n_mc = getHist<TH1D>(input_mcnumbers, "OmegaMCCrossSection/mesonCounts");

    const auto MC_Total_events = n_mc->GetEntries();
    cout << "Numer of MC input events: " << MC_Total_events << endl;


    WrapTFileInput input_data(cmd_data->getValue());
    WrapTFileInput input_mc(cmd_mc->isSet() ? cmd_mc->getValue(): cmd_data->getValue());

    // create TRint as RooFit internally creates functions/histograms, sigh...
    argc=0; // prevent TRint to parse any cmdline
    TRint app("OmegaEtaG_fit",&argc,argv,nullptr,0,true);



    unique_ptr<WrapTFileOutput> masterFile;
    if(cmd_output->isSet()) {
        // cd into masterFile upon creation
        masterFile = std_ext::make_unique<WrapTFileOutput>(cmd_output->getValue(), true);
    }


    auto lumi = getHist<TH1D>(input_lumi, "PhotonFlux/intlumicor");

    const auto Integrate= [] (const TH1* h) -> ValError {
        ValError res;
        res.v = h->IntegralAndError(1, h->GetNbinsX(), res.e);
        return res;
    };
    const auto total_lumi = Integrate(lumi);

//    const auto getCorrected = [&lumi] (WrapTFileInput& f, const string& hpath) {
//        auto h = getHist<TH2D>(f,hpath);
//        h->Divide(lumi);
//        return h->ProjectionX();
//    };

//    auto h_im_corr   = getCorrected(input_data, cmd_histpath->getValue()+datahist+"ggg_IM_taggch");
//    canvas("uncorr") << h_im_direct << endc;
//    canvas("corr")  << h_im_corr << endc;

    cosTbins_t ctbins;


    if(cmd_mode->getValue() == "global") {
        ctbins.resize(1);
        ctbins.at(0) = {NaN,NaN,{
                getHist<TH1D>(input_data, cmd_histpath->getValue()+datahist+cmd_histname->getValue()),
                getHist<TH1D>(input_mc  , cmd_histpath->getValue()+refhist+cmd_histname->getValue()),
                        n_mc->Integral(),
                                   total_lumi}};
    } else if(cmd_mode->getValue() == "cosT") {

        const auto nbins=int(n_mc->GetNbinsX());

        ctbins.resize(unsigned(nbins));
        TGraphErrors* g = new TGraphErrors(nbins);
        TGraphErrors* geff = new TGraphErrors(nbins);
        const auto SetPoint = [] (TGraphErrors& g, const int& i, const ValError& x, const ValError& y) {
            g.SetPoint(i,x.v, y.v);
            g.SetPointError(i,x.e,y.e);
        };

        for(size_t i=0;i<unsigned(nbins);++i) {
            const auto cosT = n_mc->GetBinCenter(int(i+1));
            const string basepath = std_ext::formatter() << cmd_histpath->getValue() << "/cosT_" << i;
            ctbins.at(i) = {n_mc->GetBinCenter(int(i+1)),NaN,
                    {
                        getHist<TH1D>(input_data,std_ext::formatter() << basepath << datahist << cmd_histname->getValue()),
                        getHist<TH1D>(input_mc,  std_ext::formatter() << basepath << refhist  << cmd_histname->getValue()),
                        n_mc->GetBinContent(int(1+i)),
                           total_lumi, n_mc->GetBinWidth(int(1+i))
                    }};

            const auto& fitres = ctbins.at(i);
            SetPoint(*g,    int(i), {cosT, 0.}, get<2>(fitres).vn_corr);
            SetPoint(*geff, int(i), {cosT, 0.}, get<2>(fitres).rec_eff);
        }

        auto c = new TCanvas();
        c->Divide(2,1);
        c->cd(1);
        g->Draw("AP");
        c->cd(2);
        geff->Draw("AP");

        if(cmd_ITtest->isSet()) {
            TGraph* io = new TGraph(int(ctbins.size()));
            int i=0;
            for(const auto& c : ctbins) {
                io->SetPoint(i++, get<0>(c), get<2>(c).DataToMc);
            }
            new TCanvas();
            io->Draw("AP");
            io->SetName("IOTEST");
            io->GetXaxis()->SetTitle("cos(#theta)");
            io->GetYaxis()->SetTitle("fit/input");
            gDirectory->Add(io);
        }
    }



    const auto delim = '\t';

    const auto print_cosTbins = [] (ostream& stream, const cosTbins_t& v) {
        stream << "#";
        for(const auto& h : {"cosT","E","Nsig","dNsig","Nbkg","dNbkg","RecEff","dRecEff","Ncrr","dNcorr","DataToMC","sigma","dsigma"}) {
            stream << setw(12) << h <<delim;
        }
        stream << "\n";

        for(const auto& c : v) {
            const auto& t = get<2>(c);

            stream << setw(12);
            if(isfinite(get<0>(c)))
                stream << get<0>(c);
            else stream << "";
            stream << delim;

            stream << setw(12);
            if(isfinite(get<1>(c)))
                stream << get<1>(c);
            else stream << "";
            stream << delim;

            stream << setw(12) << t.vnsig.v << delim
                 << setw(12) << t.vnsig.e << delim
                 << setw(12) << t.vnbkg.v << delim
                 << setw(12) << t.vnbkg.e << delim
                 << setw(12) << t.rec_eff.v << delim
                 << setw(12) << t.rec_eff.e << delim
                 << setw(12) << t.vn_corr.v << delim
                 << setw(12) << t.vn_corr.e << delim
                 << setw(12) << t.DataToMc << delim
                 << setw(12) << t.sigmaOmega.v << delim
                 << setw(12) << t.sigmaOmega.e
                 << "\n";
        }
        stream << endl;
    };

    print_cosTbins(cout, ctbins);

    {
        auto gnuplot = ofstream(cmd_output->getValue()+".dat");
        print_cosTbins(gnuplot, ctbins);
    }


    cout << "Total int lumi corr = " << total_lumi << endl;


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

    return 0;
}

interval<double> RoundToBin(const TH1* h, const interval<double>& i) {
    const auto lb = h->GetXaxis()->FindBin(i.Start());
    const auto hb = h->GetXaxis()->FindBin(i.Stop());
    return {h->GetXaxis()->GetBinLowEdge(lb), h->GetXaxis()->GetBinUpEdge(hb)};
}

TH1D* CutRange(const TH1* h, const interval<double>& i) {
    const auto startbin = h->GetXaxis()->FindBin(i.Start());
    const auto stopbin  = h->GetXaxis()->FindBin(i.Stop());
    TH1D* n = new TH1D("","",stopbin-startbin+1, h->GetXaxis()->GetBinLowEdge(startbin), h->GetXaxis()->GetBinUpEdge(stopbin) );

    for(int i=startbin;i<=stopbin;++i) {
        n->SetBinContent(i-startbin+1, h->GetBinContent(i));
        n->SetBinError(i-startbin+1, h->GetBinError(i));
    }
    return n;
}

FitOmegaPeak::FitOmegaPeak(const TH1 *h_data, const TH1 *h_mc, const double n_mc_input, const ValError& total_lumi, const double binwidth)
{
    auto data = h_data;
    const interval<double> fitrange = TH_ext::getBins(data->GetXaxis());
    LOG(INFO) << "Fit Range: " << fitrange;
    const auto signalregion = interval<double>::CenterWidth(ParticleTypeDatabase::Omega.Mass(), 120.0);
    LOG(INFO) << "Signal Region: " << signalregion;

    TCanvas *cfit = new TCanvas();
    cfit->SetTitle(Form("Fit: %s", h_data->GetTitle()));
    cfit->Divide(1,2);

    // define observable and ranges
    RooRealVar var_IM("IM","IM", fitrange.Start(), fitrange.Stop(), "MeV");
    var_IM.setBins(10000);
    var_IM.setRange("full",fitrange.Start(), fitrange.Stop());

    // load data to be fitted
    RooDataHist h_roo_data("h_roo_data","dataset",var_IM,h_data);

    // build shifted mc lineshape
    const double max_pos = h_data->GetBinCenter(h_data->GetMaximumBin()) - 782.0;
    RooRealVar var_IM_shift("var_IM_shift", "shift in IM", max_pos, -25.0, 25.0);
    RooProduct var_IM_shift_invert("var_IM_shift_invert","shifted IM",RooArgSet(var_IM_shift, RooConst(-1.0)));
    RooAddition var_IM_shifted("var_IM_shifted","shifted IM",RooArgSet(var_IM,var_IM_shift_invert));
    RooDataHist h_roo_mc("h_roo_mc","MC lineshape", var_IM, h_mc);
    RooHistPdf pdf_mc_lineshape("pdf_mc_lineshape","MC lineshape as PDF", var_IM_shifted, var_IM, h_roo_mc, 2);

    // build gaussian
    RooRealVar  var_gauss_sigma("gauss_sigma","width of gaussian", 9.0, 0.0, 20.0);
    RooGaussian pdf_gaussian("pdf_gaussian","Gaussian smearing", var_IM, RooConst(0.0), var_gauss_sigma);

    // build signal as convolution, note that the gaussian must be the second PDF (see documentation)
    RooFFTConvPdf pdf_signal("pdf_signal","MC_lineshape (X) gauss",var_IM, pdf_mc_lineshape, pdf_gaussian) ;

    const int polOrder = 3;
    std::vector<RooRealVar*> bkg_params; // RooRealVar cannot be copied, so create them on heap
    RooArgSet roo_bkg_params;
    for(int p=0;p<polOrder;p++) {

        bkg_params.emplace_back(new RooRealVar((
                                    "p_"+to_string(p)).c_str(), ("Bkg Par "+to_string(p)).c_str(), 0.1, -1E+10, +1E+10)); // Die Startwerte sind wichtig!!
        roo_bkg_params.add(*bkg_params.back());
    }
    RooPolynomial pdf_polbackground("Pol","Polynomial background",var_IM,roo_bkg_params);

    RooRealVar var_argus_c("argus_c","argus shape parameter",-16.0,-1000.,-0.0) ;
    RooRealVar var_argus_p("argus_p","argus shape parameter", 1.7, 0, 20) ;
    RooRealVar var_f_argus("f_argus", "argus fraction", 0.6, 0, 1);
    RooArgusBG pdf_argus("background","Argus PDF",var_IM, RooConst(990),var_argus_c, var_argus_p) ;
    RooAddPdf  pdf_background("totalbkg", "total background", RooArgList(pdf_argus, pdf_polbackground), RooArgList(var_f_argus));

    var_IM.setRange("bkg_l", fitrange.Start(), signalregion.Start());
    var_IM.setRange("bkg_r", signalregion.Stop(), fitrange.Stop());

    var_argus_p.setVal(6);
    var_argus_p.setConstant(true);

    RooRealVar nbkg = RooRealVar("nbkg","#background events", h_data->GetEntries()/2, 0, h_data->GetEntries());

    //RooExtendPdf bkg_ext = RooExtendPdf("bkg_ext","background ext",pdf_background, nbkg);
    //bkg_ext.chi2FitTo(h_roo_data, Range("bkg_l,bkg_r"), Extended());
    //   pdf_background.fitTo(h_roo_data, Range("full"), Extended());
    pdf_background.fitTo(h_roo_data, Range("bkg_l,bkg_r"), Extended(), PrintLevel(0));

    // build sum
    RooRealVar nsig = RooRealVar("nsig","#signal events", h_data->GetEntries()/2, 0, h_data->GetEntries());

    RooAddPdf pdf_sum("pdf_sum","total sum",RooArgList(pdf_signal,pdf_background),RooArgList(nsig,nbkg));

    // do the actual maximum likelihood fit
    pdf_sum.fitTo(h_roo_data,
          Extended(), Minos(RooArgSet(nsig)), SumW2Error(kTRUE), Range("full"), Save(),
          PrintLevel(0) /*,Minos(kTRUE)*/);

    // draw output, won't be shown in batch mode
    cfit->cd(1);
    RooPlot* frame = var_IM.frame();
    h_roo_data.plotOn(frame);
    //    pdf_background.plotOn(frame);
    pdf_sum.plotOn(frame, LineColor(kRed), PrintEvalErrors(-1));
    RooHist* hresid = frame->residHist();
    hresid->GetXaxis()->SetRangeUser(fitrange.Start(), fitrange.Stop());

    pdf_sum.plotOn(frame, Components(pdf_background), LineColor(kBlue), PrintEvalErrors(-1));
    pdf_sum.plotOn(frame, Components(pdf_signal), LineColor(kGreen));
    frame->Draw();

    cfit->cd(2);
    hresid->Draw();
    vnsig = nsig;
    vnbkg = nbkg;
//    rec_eff.v = h_mc->GetEntries() / n_mc_input;
//    rec_eff.e = sqrt(h_mc->GetEntries()) / n_mc_input;
    rec_eff = ValError::Statistical(h_mc->GetEntries()) / n_mc_input;

//    vn_corr.v = vnsig.v / rec_eff.v;
//    vn_corr.e = sqrt(sqr(vnsig.e / rec_eff.v) + sqr(vnsig.v/sqr(rec_eff.v)*rec_eff.e));
    vn_corr = vnsig / rec_eff;

    DataToMc = vnsig.v / h_mc->Integral();

    sigmaOmega = vn_corr / total_lumi / BRomegaPi0ggg /binwidth;

   // fr->Print();
}

ostream& operator<<(ostream &s, const FitOmegaPeak &f)
{
    s << "[NSig=" << f.vnsig
      << " Nbkg=" << f.vnbkg
      << " Npar=" << f.numParams
      << " chi2dof=" << f.chi2ndf
      << " RecEff=" << f.rec_eff
      << " N_corr=" << f.vn_corr
      << " Data/MC=" << f.DataToMc
      << "]";
    return s;
}
