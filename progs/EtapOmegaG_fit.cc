#include "base/Logger.h"

#include "tclap/CmdLine.h"
#include "tclap/ValuesConstraintExtra.h"
#include "base/interval.h"
#include "base/WrapTFile.h"
#include "base/std_ext/string.h"
#include "base/std_ext/system.h"
#include "base/std_ext/memory.h"
#include "base/std_ext/math.h"
#include "base/ParticleType.h"

#include "analysis/plot/RootDraw.h"
#include "analysis/utils/ParticleTools.h"
#include "root-addons/analysis_codes/Math.h"
#include "expconfig/ExpConfig.h"
#include "base/Detector_t.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TLatex.h"

#include "TSystem.h"
#include "TRint.h"

#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooArgusBG.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooHistPdf.h"
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
#include "RooHist.h"
#include "RooPlotable.h"

#include "APLCON.hpp"

auto debug = false;

using namespace ant;
using namespace std;
using namespace RooFit;

struct fit_params_t {
    interval<double> signal_region{920, 990};
    int nSamplingBins{10000};
    int interpOrder{2};
    double ymax{160};

    unsigned TaggCh = 0;
    double Eg = std_ext::NaN;

    TH1D* h_mc = nullptr;
    TH1D* h_data = nullptr;

};

struct fit_return_t : ant::root_drawable_traits {

    fit_params_t p;

    RooFitResult* fitresult = nullptr;
    double chi2ndf = std_ext::NaN;
    double peakpos = std_ext::NaN;

    int numParams() {
        return fitresult->floatParsFinal().getSize();
    }

    double residualSignalIntegral() {
        auto& h = residual;
        return h->Integral(h->GetXaxis()->FindBin(p.signal_region.Start()),
                           h->GetXaxis()->FindBin(p.signal_region.Stop()))
                /h->getNominalBinWidth();
    }

    const RooRealVar& getNfit() const {
        auto& pars = fitresult->floatParsFinal();
        return dynamic_cast<const RooRealVar&>(*pars.at(pars.index("nsig")));
    }

    RooPlot* fitplot = nullptr;
    TPaveText* lbl = nullptr;
    RooHist* residual = nullptr;

    virtual void Draw(const string& option) const override
    {
        auto& nsig = getNfit();
        auto& pars = fitresult->floatParsFinal();
        auto& sigma = dynamic_cast<const RooRealVar&>(*pars.at(pars.index("sigma")));
        auto& shift = dynamic_cast<const RooRealVar&>(*pars.at(pars.index("var_IM_shift")));

        (void)option;

        auto lbl = new TPaveText();
        lbl->SetX1NDC(0.65);
        lbl->SetX2NDC(0.98);
        lbl->SetY1NDC(0.5);
        lbl->SetY2NDC(0.95);
        lbl->SetBorderSize(0);
        lbl->SetFillColor(kWhite);
        lbl->AddText(static_cast<string>(std_ext::formatter() << setprecision(1) << fixed << "E_{#gamma} = " << p.Eg << " MeV").c_str());
        lbl->AddText(static_cast<string>(std_ext::formatter() << setprecision(0) << fixed << "N_{sig} = " << nsig.getVal() << " #pm " << nsig.getError()).c_str());
        lbl->AddText(static_cast<string>(std_ext::formatter() << setprecision(2) << fixed << "#chi^{2}_{red} = " << chi2ndf).c_str());
        lbl->AddText(static_cast<string>(std_ext::formatter() << setprecision(1) << fixed << "#sigma = " << sigma.getVal() << " MeV").c_str());
        lbl->AddText(static_cast<string>(std_ext::formatter() << setprecision(1) << fixed << "#delta = " << shift.getVal() << " MeV").c_str());

        if(debug) {
            std_ext::formatter extra;
            extra << "Status: ";
            for(unsigned i=0;i<fitresult->numStatusHistory();i++)
                extra << fitresult->statusCodeHistory(i);
            extra << " TaggCh=" << p.TaggCh;
            lbl->AddText(static_cast<string>(extra).c_str());

            if(fitresult->status())
                lbl->SetTextColor(kRed);

            if(nsig.getError()>200)
                lbl->SetTextColor(kRed);

            if(chi2ndf > 2)
                lbl->SetTextColor(kRed);
        }

        fitplot->SetMinimum(0);
        fitplot->SetMaximum(p.ymax);
        fitplot->SetTitle("");
        fitplot->Draw();
        lbl->Draw();
        gPad->SetTopMargin(0.01);
        gPad->SetRightMargin(0.003);
        gPad->SetLeftMargin(0.07);
    }

    friend ostream& operator<<(ostream& s, const fit_return_t& o) {
        const auto options = "v";
        o.fitresult->printStream(s,o.fitresult->defaultPrintContents(options),o.fitresult->defaultPrintStyle(options));
        return s;
    }
};

fit_return_t doFit_Ref(const fit_params_t& p) {

    fit_return_t r;
    r.p = p; // remember params

    const auto calcIMThresh = [] (double Eg) {
        const auto mp = ParticleTypeDatabase::Proton.Mass();
        return std::sqrt(std_ext::sqr(mp) + 2*mp*Eg) - mp;
    };
    const auto threshold = calcIMThresh(p.Eg);

    // define observable and ranges
    RooRealVar var_IM("IM","IM", p.h_data->GetXaxis()->GetXmin(), p.h_data->GetXaxis()->GetXmax(), "MeV");
    var_IM.setBins(p.nSamplingBins);
    var_IM.setRange("full",var_IM.getMin(),threshold+2); // for close-to-threshold tagger bins

    // load data to be fitted
    RooDataHist h_roo_data("h_roo_data","dataset",var_IM,p.h_data);

    // build shifted mc lineshape
    RooRealVar var_IM_shift("var_IM_shift", "shift in IM", 0.0, -10.0, 10.0);
    RooProduct var_IM_shift_invert("var_IM_shift_invert","shifted IM",RooArgSet(var_IM_shift, RooConst(-1.0)));
    RooAddition var_IM_shifted("var_IM_shifted","shifted IM",RooArgSet(var_IM,var_IM_shift_invert));
    RooDataHist h_roo_mc("h_roo_mc","MC lineshape", var_IM, p.h_mc);
    RooHistPdf pdf_mc_lineshape("pdf_mc_lineshape","MC lineshape as PDF", var_IM_shifted, var_IM, h_roo_mc, p.interpOrder);

    // build detector resolution smearing

    RooRealVar  var_sigma("sigma","detector resolution",  2.0, 0.0, 10.0);
    RooGaussian pdf_smearing("pdf_smearing","Single Gaussian", var_IM, RooConst(0.0), var_sigma);

    // build signal as convolution, note that the gaussian must be the second PDF (see documentation)
    RooFFTConvPdf pdf_signal("pdf_signal","MC_lineshape (X) gauss", var_IM, pdf_mc_lineshape, pdf_smearing);

    // build background (chebychev or argus?)

    //    const int polOrder = 6;
    //    std::vector<std::unique_ptr<RooRealVar>> bkg_params; // RooRealVar cannot be copied, so create them on heap
    //    RooArgSet roo_bkg_params;
    //    for(int p=0;p<polOrder;p++) {
    //        bkg_params.emplace_back(std_ext::make_unique<RooRealVar>((
    //                                    "p_"+to_string(p)).c_str(), ("Bkg Par "+to_string(p)).c_str(), 0.0, -10.0, 10.0)
    //                                );
    //        roo_bkg_params.add(*bkg_params.back());
    //    }
    //    RooChebychev pdf_background("chebychev","Polynomial background",var_IM,roo_bkg_params);
    //    var_IM.setRange("bkg_l", var_IM.getMin(), 930);
    //    var_IM.setRange("bkg_r", 990, 1000);
    //    pdf_background.fitTo(h_roo_data, Range("bkg_l,bkg_r"), Extended()); // using Range(..., ...) does not work here (bug in RooFit, sigh)



    RooRealVar argus_cutoff("argus_cutoff","argus pos param", threshold);
    RooRealVar argus_shape("argus_shape","argus shape param", -5, -25.0, 0.0);
    RooRealVar argus_p("argus_p","argus p param", 0.5);
    RooArgusBG pdf_background("argus","bkg argus",var_IM,argus_cutoff,argus_shape,argus_p);

    // build sum
    RooRealVar nsig("nsig","#signal events", 3e3, 0, 1e6);
    RooRealVar nbkg("nbkg","#background events", 3e3, 0, 1e6);
    RooAddPdf pdf_sum("pdf_sum","total sum",RooArgList(pdf_signal,pdf_background),RooArgList(nsig,nbkg));

    // do some pre-fitting to obtain better starting values, make sure function is non-zero in range
//    var_IM.setRange("nonzero",var_IM.getMin(), threshold-5);
//    pdf_sum.chi2FitTo(h_roo_data, Range("nonzero"), PrintLevel(-1)); // using Range(..., ...) does not work here (bug in RooFit, sigh)


    // do the actual maximum likelihood fit
    // use , Optimize(false), Strategy(2) for double gaussian...?!
    r.fitresult = pdf_sum.fitTo(h_roo_data, Extended(), SumW2Error(kTRUE), Range("full"), Save(), PrintLevel(-1));

    // draw output and remember pointer
    r.fitplot = var_IM.frame();

    h_roo_data.plotOn(r.fitplot);
    //    pdf_sum.plotOn(frame, LineColor(kRed), VisualizeError(*fr));

    // need to figure out chi2nds and stuff after plotting data and finally fitted pdf_sum
    // also the residHist must be created here (and rememebered for later use)
    pdf_sum.plotOn(r.fitplot, LineColor(kRed));
    r.chi2ndf = r.fitplot->chiSquare(r.numParams());
    auto pdf_sum_tf = pdf_sum.asTF(var_IM);
    r.peakpos = pdf_sum_tf->GetMaximumX(p.signal_region.Start(), p.signal_region.Stop());
    r.residual = r.fitplot->residHist();

    pdf_sum.plotOn(r.fitplot, Components(pdf_background), LineColor(kBlue));
    pdf_sum.plotOn(r.fitplot, Components(pdf_signal), LineColor(kGreen));

    return r;
}

// use APLCON to calculate the total sum with error propagation
// Value, Sigma, Pull
struct N_t {
    N_t(const RooRealVar& var) : Value(var.getVal()), Sigma(var.getError()) {}
    N_t(double v=0, double s=0) : Value(v), Sigma(s) {}

    double Value;
    double Sigma;
    double Pull = std_ext::NaN;

    template<std::size_t N>
    std::tuple<double&> linkFitter() noexcept {
        // the following get<N> assumes this order of indices
        static_assert(APLCON::ValueIdx==0,"");
        static_assert(APLCON::SigmaIdx==1,"");
        static_assert(APLCON::PullIdx ==2,"");
        // the extra std::tie around std::get is for older compilers...
        return std::tie(std::get<N>(std::tie(Value, Sigma, Pull)));
    }

    friend ostream& operator<<(ostream& s, const N_t& o) {
        return s << o.Value << "+/-" << o.Sigma << "(" << 100.0*o.Sigma/o.Value << "%)";
    }
};

int main(int argc, char** argv) {

    SetupLogger();

    TCLAP::CmdLine cmd("EtapOmegaG_fit", ' ', "0.1");
    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"int");
    auto cmd_batchmode = cmd.add<TCLAP::MultiSwitchArg>("b","batch","Run in batch mode (no ROOT shell afterwards)",false);
    auto cmd_debug = cmd.add<TCLAP::MultiSwitchArg>("","debug","Enable debug mode",false);
    auto cmd_output = cmd.add<TCLAP::ValueArg<string>>("o","output","Output file",false,"","filename");
    TCLAP::ValuesConstraintExtra<decltype(ExpConfig::Setup::GetNames())> allowedsetupnames(ExpConfig::Setup::GetNames());
    auto cmd_setup  = cmd.add<TCLAP::ValueArg<string>>("s","setup","Choose setup by name",true,"", &allowedsetupnames);

    auto cmd_data = cmd.add<TCLAP::ValueArg<string>>("","data","Data input",true,"","rootfile");
    auto cmd_mc = cmd.add<TCLAP::ValueArg<string>>("","mc","MC signal/reference input",true,"","rootfile");

    cmd.parse(argc, argv);
    if(cmd_verbose->isSet()) {
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());
    }
    debug = cmd_debug->isSet();

    RooMsgService::instance().setGlobalKillBelow(debug ? RooFit::WARNING : RooFit::ERROR);

    ExpConfig::Setup::SetByName(cmd_setup->getValue());
    auto Tagger = ExpConfig::Setup::GetDetector<TaggerDetector_t>();

    const string ref_prefix   = "EtapOmegaG_plot_Ref";
    const string ref_histpath = ref_prefix+"/DiscardedEk=0/KinFitProb>0.02";
    const string ref_histname = "h_IM_2g_TaggCh";

    TH2D* ref_data;
    TH2D* ref_mc;
    TH1D* ref_mctrue_generated;

    WrapTFileInput input_data(cmd_data->getValue()); // keep it open
    {
        const string histpath = ref_histpath+"/h/Data/"+ref_histname;
        if(!input_data.GetObject(histpath, ref_data)) {
            LOG(ERROR) << "Cannot find " << histpath;
            return EXIT_FAILURE;
        }
    }

    WrapTFileInput input_mc(cmd_mc->getValue()); // keep it open
    {
        const string histpath = ref_histpath+"/h/Ref/"+ref_histname;
        if(!input_mc.GetObject(histpath, ref_mc)) {
            LOG(ERROR) << "Cannot find " << histpath;
            return EXIT_FAILURE;
        }
    }
    {
        const string histpath = ref_prefix+"/h_mctrue_generated";
        if(!input_mc.GetObject(histpath, ref_mctrue_generated)) {
            LOG(ERROR) << "Cannot find " << histpath;
            return EXIT_FAILURE;
        }
    }

    // create TRint as RooFit internally creates functions/histograms, sigh...
    argc=0; // prevent TRint to parse any cmdline
    TRint app("EtapOmegaG_plot",&argc,argv,nullptr,0,true);

    unique_ptr<WrapTFileOutput> masterFile;
    if(cmd_output->isSet()) {
        // cd into masterFile upon creation
        masterFile = std_ext::make_unique<WrapTFileOutput>(cmd_output->getValue(), true);
    }


    ant::canvas("EtapOmegaG_fit: Input")
            << drawoption("colz")
            << ref_mc << ref_data << ref_mctrue_generated
            << endc;


    ant::canvas c_plots_data("EtapOmegaG_fit: Plots Data");

    std::vector<N_t> Ns_fit;
    std::vector<N_t> Ns_effcorr;

    const auto maxTaggCh = 1; // should be 40 or 39

    for(int taggch=maxTaggCh;taggch>=0;taggch--) {
        LOG(INFO) << "Fitting TaggCh=" << taggch;
        fit_params_t p;
        p.TaggCh = taggch;
        p.Eg = Tagger->GetPhotonEnergy(taggch);

        // fit MC lineshape to data
        const auto taggbin = taggch+1;
        p.h_mc   = ref_mc->ProjectionX("h_mc",taggbin,taggbin);
        p.h_data = ref_data->ProjectionX("h_data",taggbin,taggbin);
        auto r_data = doFit_Ref(p);

        // determine efficiency corrected N_effcorr = N_fit/efficiency = N_fit * mc_generated/mc_reco;
        N_t N_fit(r_data.getNfit());
        Ns_fit.emplace_back(N_fit); // store here as APLCON is going to change it
        N_t N_effcorr;
        {
            N_t N_mcreco;
            N_mcreco.Value = p.h_mc->IntegralAndError(1, p.h_mc->GetNbinsX(), N_mcreco.Sigma, ""); // take binwidth into account?
            N_t N_mcgen(ref_mctrue_generated->GetBinContent(taggbin), ref_mctrue_generated->GetBinError(taggbin));

            APLCON::Fitter<N_t, N_t, N_t, N_t> calcN_effcorr;
            calcN_effcorr.DoFit(N_fit, N_mcreco, N_mcgen, N_effcorr,
                                [] (const N_t& N_fit, const N_t& N_mcreco, const N_t& N_mcgen, const N_t& N_effcorr) {
                return N_effcorr.Value - N_fit.Value * N_mcgen.Value / N_mcreco.Value;
            });

            if(debug) {
                LOG(INFO) << "TaggCh=" << taggch
                          << " N_fit=" << Ns_fit.back()
                          << " N_effcorr=" << N_effcorr;
            }
        }

        // save/plot values
        Ns_effcorr.emplace_back(N_effcorr);

        c_plots_data << r_data;
        if(debug) {
            LOG(INFO) << r_data;
        }
    }
    c_plots_data << endc;

    // sum up the N_data and N_effcorr
    N_t Nsum; // sigma=0 means unmeasured
    LOG(INFO) << Ns_fit;
    {
        APLCON::Fitter<std::vector<N_t>, N_t> sum_Nsig_data;
        sum_Nsig_data.DoFit(Ns_fit, Nsum, [] (const vector<N_t>& N, const N_t& Nsum) {
            double sum = 0.0;
            for(auto& n : N)
                sum += n.Value;
            return Nsum.Value - sum;
        });
    }

    LOG(INFO) << "THE TOTAL SUM IS: " << Nsum;

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
