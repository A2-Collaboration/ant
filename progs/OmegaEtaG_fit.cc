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
#include "expconfig/setups/Setup.h"
#include "expconfig/setups/SetupRegistry.h"
#include "expconfig/detectors/Tagger.h"

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
#include "base/PhysicsMath.h"
#include "TPavesText.h"

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
    int ndf = -1;
    double chi2ndf = std_ext::NaN;

    ValError rec_eff = std_ext::NaN;
    ValError vn_corr = std_ext::NaN;
    ValError sigmaOmega = std_ext::NaN;
    ValError sigma = std_ext::NaN;
    ValError argus_c = std_ext::NaN;
    ValError argus_p = std_ext::NaN;
    ValError argus_f = std_ext::NaN;
    ValError mshift = std_ext::NaN;

    double nMC = std_ext::NaN;
    double nMCInput = std_ext::NaN;

    FitOmegaPeak() {}
    FitOmegaPeak(const TH1* data, const TH1* mc_shape, const double n_mc_input, const ValError& total_lumi, const double binwidth=1.0, const interval<double> fitrange={670,900}, TVirtualPad* pad=nullptr, const string &title="");
    FitOmegaPeak(const FitOmegaPeak&) = default;
    FitOmegaPeak(FitOmegaPeak&&) = default;
    FitOmegaPeak& operator=(const FitOmegaPeak&) = default;
    FitOmegaPeak& operator=(FitOmegaPeak&&) = default;

    friend ostream& operator<<(ostream& s, const FitOmegaPeak& argus_f);
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

using cosTbins_t = vector<tuple<double,interval<double>,FitOmegaPeak>>;

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
    auto cmd_cosTslice = cmd.add<TCLAP::ValueArg<int>>("", "cosTslice","",false,-1,"slice");
    auto cmd_Taggslice = cmd.add<TCLAP::ValueArg<int>>("", "Taggslice","",false,-1,"slice");

  //  auto cmd_range     = cmd.add<TCLAP::ValueArg<interval<double>>>("","range","Fit range",false,"","range");

    cmd.parse(argc, argv);
    if(cmd_verbose->isSet()) {
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());
    }

    const interval<double> fitrange = {670.0,900.0};

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

    const auto Integrate= [] (const TH1* h, int start=1, int end=-1) -> ValError {
        if(end==-1)
            end=h->GetNbinsX();
        ValError res;
        res.v = h->IntegralAndError(start, end, res.e);
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
        ctbins.reserve(1);
        auto h_data =  getHist<TH1D>(input_data, cmd_histpath->getValue()+datahist+cmd_histname->getValue());
        ctbins.push_back({NaN,interval<double>(NaN,NaN),{
               h_data,
                getHist<TH1D>(input_mc  , cmd_histpath->getValue()+refhist+cmd_histname->getValue()),
                        n_mc->Integral(),
                                   total_lumi,1.0,TH_ext::getBins(h_data->GetXaxis())}});
    } else if(cmd_mode->getValue() == "cosT") {

        const auto nbins=int(n_mc->GetNbinsX());

        ctbins.reserve(unsigned(nbins));
        TGraphErrors* g = new TGraphErrors(nbins);
        TGraphErrors* geff = new TGraphErrors(nbins);
        const auto SetPoint = [] (TGraphErrors& g, const int& i, const ValError& x, const ValError& y) {
            g.SetPoint(i,x.v, y.v);
            g.SetPointError(i,x.e,y.e);
        };

        for(size_t i=0;i<unsigned(nbins);++i) {
            const auto cosT = n_mc->GetBinCenter(int(i+1));
            const string basepath = std_ext::formatter() << cmd_histpath->getValue() << "/cosT_" << i;
            const auto h_data = getHist<TH1D>(input_data,std_ext::formatter() << basepath << datahist << cmd_histname->getValue());
            const auto h_mc = getHist<TH1D>(input_mc,  std_ext::formatter() << basepath << refhist  << cmd_histname->getValue());
            ctbins.push_back({n_mc->GetBinCenter(int(i+1)),{NaN,NaN},
                    {
                        h_data,
                        h_mc,
                        n_mc->GetBinContent(int(1+i)),
                           total_lumi, n_mc->GetBinWidth(int(1+i)),TH_ext::getBins(h_data->GetXaxis())
                    }});

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

//        if(cmd_ITtest->isSet()) {
//            TGraph* io = new TGraph(int(ctbins.size()));
//            int i=0;
//            for(const auto& c : ctbins) {
//                io->SetPoint(i++, get<0>(c), get<2>(c).nMC);
//            }
//            new TCanvas();
//            io->Draw("AP");
//            io->SetName("IOTEST");
//            io->GetXaxis()->SetTitle("cos(#theta)");
//            io->GetYaxis()->SetTitle("fit/input");
//            gDirectory->Add(io);
//        }
    } else if(cmd_mode->getValue() == "cosTE") {

        ExpConfig::Setup::SetByName("Setup_2014_10_EPT_Prod");
        auto tagger = ExpConfig::Setup::GetDetector<TaggerDetector_t>();

        constexpr int ntaggergroup = 4;

        const auto EWindow = [&tagger, &ntaggergroup] (const int tagger_group) {
            const auto clip = [&tagger] (const unsigned t) { return t < tagger->GetNChannels() ? t :tagger->GetNChannels() -1; };
            const auto t1 = unsigned(tagger_group*ntaggergroup);
            const auto t2 = unsigned(clip((tagger_group+1)*ntaggergroup-1));

            return interval<double>(
                tagger->GetPhotonEnergy(t1) - tagger->GetPhotonEnergyWidth(t1)/2.0,
                tagger->GetPhotonEnergy(t2) + tagger->GetPhotonEnergyWidth(t2)/2.0
            );

        };

//        const auto getTagSlice = [&ntaggergroup] (const TH2D* h, int n) {
//            return h->ProjectionX(Form("%s_%d",h->GetName(), n),n,n+ntaggergroup-1);
//        };

        auto n_mc_input_2d = getHist<TH2D>(input_mcnumbers, "OmegaMCCrossSection/mesonCounts_taggch");



        const auto nc = n_mc_input_2d->GetNbinsX();
        const auto n_tagger_bins = n_mc_input_2d->GetNbinsY();
        const auto n_tagger_groups = int(ceil(n_tagger_bins / double(ntaggergroup)));



        //canvas->Divide(nc,n_tagger_bins);

        ctbins.reserve( unsigned(nc * n_tagger_groups ));

        const auto cosTslicelimit = cmd_cosTslice->getValue();
        const auto Taggslicelimit = cmd_Taggslice->getValue();
        TCanvas* canvas = nullptr;

        for(int c=1;c<=nc; ++c) {

            if(cosTslicelimit >= 0 && (c-1)!= cosTslicelimit)
                continue;

            const auto cosT = n_mc_input_2d->GetXaxis()->GetBinCenter(int(c));
            const auto cosT_binwidth = n_mc_input_2d->GetXaxis()->GetBinWidth(int(c));

            const string basepath = std_ext::formatter() << cmd_histpath->getValue() << "/cosT_" << c-1;

            const auto h_data2d = getHist<TH2D>(input_data, std_ext::formatter() << basepath << datahist << cmd_histname->getValue());
            const auto h_mc2d   = getHist<TH2D>(input_mc,   std_ext::formatter() << basepath << refhist  << cmd_histname->getValue());
            const auto n_mc_input = n_mc_input_2d->ProjectionY("",c,c);

            for(int tagger_group=0; tagger_group<n_tagger_groups; tagger_group++) {

                if(Taggslicelimit >= 0 && tagger_group!= Taggslicelimit)
                    continue;

                const auto tagger_bin = tagger_group * ntaggergroup;

                auto h_data_slice = h_data2d->ProjectionX("data_slice",tagger_bin,tagger_bin+ntaggergroup-1);
                auto h_mc_slice   = h_mc2d->ProjectionX(  "mc_slice",  tagger_bin,tagger_bin+ntaggergroup-1);

                h_data_slice->GetXaxis()->SetRangeUser(fitrange.Start(), fitrange.Stop());
                h_mc_slice->GetXaxis()->SetRangeUser(fitrange.Start(), fitrange.Stop());

                const auto Eg = EWindow(tagger_group);

                const auto lumi_slice = Integrate(lumi,tagger_bin,tagger_bin+ntaggergroup-1);

                //auto pad = canvas->cd(1 + tagger_group + (c-1)*n_tagger_bins);
                const string title = formatter() << "W=" << round(math::W(Eg.Center(), ParticleTypeDatabase::Proton)*100.0)/100.0 << " MeV cos(#theta)_{cm}=" << cosT;

                delete canvas;
                canvas = new TCanvas();
                canvas->SetCanvasSize(600,600);

                ctbins.push_back({cosT,Eg,
                                     {
                                         h_data_slice,
                                         h_mc_slice,
                                         n_mc_input->Integral(int(tagger_bin),int(tagger_bin)+ntaggergroup-1),
                                         lumi_slice,
                                         cosT_binwidth,
                                         TH_ext::getBins(h_data_slice->GetXaxis()),
                                         canvas,
                                        title
                                     }});

                const string fname = formatter() << "W=" << round(math::W(Eg.Center(), ParticleTypeDatabase::Proton)*100.0)/100.0 << "cosT=" << cosT;
                canvas->SaveMultiImages(fname.c_str());
            }


        }

    } else {
        LOG(FATAL) << "invalid mode: " << cmd_mode->getValue();
    }


    const auto delim = '\t';

    const auto print_cosTbins = [] (ostream& stream, const cosTbins_t& v) {
        stream << "#";
        for(const auto& h : {"cosT","Emin","Emax", "Ecenter", "Nsig","dNsig","Nbkg","dNbkg",
            "RecEff","dRecEff","Ncrr","dNcorr","nMC","nMCInput","sigma","dsigma",
            "chi2dof","mshift","gauss","c","p","f"
    }) {
            stream << setw(12) << h <<delim;
        }
        stream << "\n";

        for (const auto& c : v) {
            const auto& t = get<2>(c);

            stream << setw(12) << get<0>(c) << delim;

            const auto& Eg = get<1>(c);

            stream << setw(12) << Eg.Start()   << delim;
            stream << setw(12) << Eg.Stop()    << delim;
            stream << setw(12) << Eg.Center()  << delim;

            stream
                 << setw(12) << t.vnsig.v      << delim
                 << setw(12) << t.vnsig.e      << delim
                 << setw(12) << t.vnbkg.v      << delim
                 << setw(12) << t.vnbkg.e      << delim
                 << setw(12) << t.rec_eff.v    << delim
                 << setw(12) << t.rec_eff.e    << delim
                 << setw(12) << t.vn_corr.v    << delim
                 << setw(12) << t.vn_corr.e    << delim
                 << setw(12) << t.nMC          << delim
                 << setw(12) << t.nMCInput     << delim
                 << setw(12) << t.sigmaOmega.v << delim
                 << setw(12) << t.sigmaOmega.e << delim
                 << setw(12) << t.chi2ndf      << delim
                 << setw(12) << t.mshift.v     << delim
                 << setw(12) << t.sigma.v      << delim
                 << setw(12) << t.argus_c.v    << delim
                 << setw(12) << t.argus_p.v    << delim
                 << setw(12) << t.argus_f.v    << delim
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

FitOmegaPeak::FitOmegaPeak(const TH1 *h_data, const TH1 *h_mc, const double n_mc_input, const ValError& total_lumi, const double binwidth, const interval<double> fitrange, TVirtualPad *pad, const string& title)
{

    LOG(INFO) << "Fit Range: " << fitrange;
    const auto signalregion = interval<double>::CenterWidth(ParticleTypeDatabase::Omega.Mass(), 120.0);
    LOG(INFO) << "Signal Region: " << signalregion;


    if(!pad) {
        pad = new TCanvas();
        pad->SetCanvasSize(600,600);
    }

    pad->SetTitle(Form("Fit: %s", h_data->GetTitle()));
    pad->Divide(1,2);

    {
        auto p = pad->cd(1);
        p->SetMargin(0.10f,0.05f,0.0f,0.10f);
        p->SetPad(0.0,0.3,1.0,1.0);
        p->SetTicks(1,1);
    }
    {
        auto p = pad->cd(2);
        p->SetMargin(0.10f,0.05f,0.25f,0.00f);
        p->SetPad(0.0,0.0,1.0,0.3);
        p->SetTicks(1,1);
    }
    pad->cd(1);

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

    const auto NTotal = h_data->Integral();

    RooRealVar nbkg = RooRealVar("nbkg","#background events", NTotal/2, 0, 2*NTotal);

    //RooExtendPdf bkg_ext = RooExtendPdf("bkg_ext","background ext",pdf_background, nbkg);
    //bkg_ext.chi2FitTo(h_roo_data, Range("bkg_l,bkg_r"), Extended());
    //   pdf_background.fitTo(h_roo_data, Range("full"), Extended());
    pdf_background.fitTo(h_roo_data, Range("bkg_l,bkg_r"), PrintLevel(0));

    // build sum
    RooRealVar nsig = RooRealVar("nsig","#signal events", NTotal/2, 0, 2*NTotal);

    RooAddPdf pdf_sum("pdf_sum","total sum",RooArgList(pdf_signal,pdf_background),RooArgList(nsig,nbkg));

    // do the actual maximum likelihood fit
    const auto res = pdf_sum.fitTo(h_roo_data,
          Extended(), Minos(RooArgSet(nsig)), SumW2Error(kTRUE), Range("full"), Save(),
          PrintLevel(0) /*,Minos(kTRUE)*/);

    // draw output, won't be shown in batch mode
    pad->cd(1);
    RooPlot* frame = var_IM.frame();
    h_roo_data.plotOn(frame);
    frame->GetYaxis()->SetLabelSize(0.05f);
    frame->GetYaxis()->SetTitleSize(0.05f);
    frame->GetXaxis()->SetNdivisions(504);
    frame->GetYaxis()->SetNdivisions(504);
    frame->GetXaxis()->SetRangeUser(fitrange.Start(), fitrange.Stop());
    frame->SetTitle(title.c_str());

    auto p = new TPaveText();
    p->SetX1NDC(0.6);
    p->SetX2NDC(0.98);
    p->SetY1NDC(0.45);
    p->SetY2NDC(0.9);
    p->SetTextSize(0.04f);

    const auto addLine = [] (TPaveText& p, const RooRealVar& v, const string& name="") {
        p.InsertText(Form("%s = %.3f#pm%.3f", name.empty() ? v.GetName() : name.c_str(), v.getValV(), v.getError()));
    };


    //    pdf_background.plotOn(frame);
    pdf_sum.plotOn(frame, LineColor(kRed), PrintEvalErrors(-1));
    RooHist* hresid = frame->residHist();
    hresid->SetTitle("");
    hresid->GetXaxis()->SetRangeUser(fitrange.Start(), fitrange.Stop());
    hresid->GetXaxis()->SetTitle("m(#pi^{0}#gamma) [MeV]");
    hresid->GetXaxis()->SetLabelSize(0.12f);
    hresid->GetXaxis()->SetTitleSize(0.12f);
    hresid->GetXaxis()->SetNdivisions(504);
    hresid->GetXaxis()->SetTickLength(.1f);
    hresid->GetYaxis()->SetNdivisions(404);
    hresid->GetYaxis()->SetLabelSize(0.12f);

    pdf_sum.plotOn(frame, Components(pdf_background), LineColor(kBlue), PrintEvalErrors(-1));
    pdf_sum.plotOn(frame, Components(pdf_signal), LineColor(kGreen));
    frame->Draw();
    pdf_sum.paramOn(frame);
    chi2ndf = frame->chiSquare(res->floatParsFinal().getSize());

    p->InsertText(Form("#chi^{2}/dof = %.3f", chi2ndf));
    addLine(*p, var_IM_shift,    "m - m_{#omega}");
    addLine(*p, var_gauss_sigma, "#sigma");
    addLine(*p, var_argus_c,     "c");
    addLine(*p, var_argus_p,     "p");
    addLine(*p, var_f_argus,     "f");
    addLine(*p, nsig,            "n_{sig}");
    addLine(*p, nbkg,            "n_{bkg}");
    p->Draw();

    pad->cd(2);
    hresid->Draw();

    vnsig = nsig;
    vnbkg = nbkg;

    nMCInput =  n_mc_input;
    nMC = h_mc->Integral();
    rec_eff = ValError::Statistical(nMC) / n_mc_input;

    vn_corr = vnsig / rec_eff;

    sigma = var_gauss_sigma;
    argus_c = var_argus_c;
    this->argus_p = var_argus_p;
    argus_f = var_f_argus;
    mshift = var_IM_shift;



    sigmaOmega = vn_corr / total_lumi / BRomegaPi0ggg /binwidth;

}

ostream& operator<<(ostream &s, const FitOmegaPeak &f)
{
    s << "[NSig="    << f.vnsig
      << " Nbkg="    << f.vnbkg
      << " Npar="    << f.numParams
      << " chi2dof=" << f.chi2ndf
      << " RecEff="  << f.rec_eff
      << " N_corr="  << f.vn_corr
      << " nMC="     << f.nMC
      << " nMCInput=" << f.nMCInput
      << "]";
    return s;
}
