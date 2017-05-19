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
#include "TGraphErrors.h"
#include "TMultiGraph.h"

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

// helper structs to pass parameters for fitting around
struct fit_params_t {
    interval<double> signal_region{920, 990};
    int nSamplingBins{10000};
    int interpOrder{4};
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
    double threshold = std_ext::NaN;

    int numParams() {
        return fitresult->floatParsFinal().getSize();
    }

    double residualSignalIntegral() {
        auto& h = residual;
        return h->Integral(h->GetXaxis()->FindBin(p.signal_region.Start()),
                           h->GetXaxis()->FindBin(p.signal_region.Stop()))
                /h->getNominalBinWidth();
    }

    N_t getNfit() const {
        auto& pars = fitresult->floatParsFinal();
        return dynamic_cast<const RooRealVar&>(*pars.at(pars.index("nsig")));
    }

    RooPlot* fitplot = nullptr;
    RooHist*  h_data = nullptr;
    RooCurve* f_sum = nullptr;
    RooCurve* f_sig = nullptr;
    RooCurve* f_bkg = nullptr;

    RooHist* residual = nullptr;
    N_t N_effcorr;

    virtual void Draw(const string& option) const override
    {
        auto nsig = getNfit();
        auto& pars = fitresult->floatParsFinal();
        auto& sigma = dynamic_cast<const RooRealVar&>(*pars.at(pars.index("sigma")));
        auto& shift = dynamic_cast<const RooRealVar&>(*pars.at(pars.index("x_shift")));

        (void)option;

        auto lbl = new TPaveText();
        lbl->SetX1NDC(0.65);
        lbl->SetX2NDC(0.98);
        lbl->SetY1NDC(0.5);
        lbl->SetY2NDC(0.95);
        lbl->SetBorderSize(0);
        lbl->SetFillColor(kWhite);
        lbl->AddText(static_cast<string>(std_ext::formatter() << setprecision(1) << fixed << "E_{#gamma} = " << p.Eg << " MeV").c_str());
        lbl->AddText(static_cast<string>(std_ext::formatter() << setprecision(0) << fixed << "N = " << nsig.Value << " #pm " << nsig.Sigma).c_str());
        lbl->AddText(static_cast<string>(std_ext::formatter() << setprecision(0) << fixed << "N/#varepsilon = " << N_effcorr.Value << " #pm " << N_effcorr.Sigma).c_str());
        lbl->AddText(static_cast<string>(std_ext::formatter() << setprecision(2) << fixed << "#chi^{2}_{red} = " << chi2ndf).c_str());
        lbl->AddText(static_cast<string>(std_ext::formatter() << setprecision(1) << fixed << "#sigma = " << sigma.getVal() << " MeV").c_str());
        lbl->AddText(static_cast<string>(std_ext::formatter() << setprecision(1) << fixed << "#delta = " << shift.getVal() << " MeV").c_str());

        if(debug) {
            std_ext::formatter extra;
            extra << "Status: ";
            for(unsigned i=0;i<fitresult->numStatusHistory();i++) {
                auto code = fitresult->statusCodeHistory(i);
                if(code != 0)
                    lbl->SetTextColor(kRed);
                extra << code;
            }
            extra << " TaggCh=" << p.TaggCh;
            lbl->AddText(static_cast<string>(extra).c_str());

            if(fitresult->status())
                lbl->SetTextColor(kRed);

            if(nsig.Sigma>200)
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
    r.threshold = calcIMThresh(p.Eg);

    // define observable and ranges
    RooRealVar x("IM","IM", p.h_data->GetXaxis()->GetXmin(), p.h_data->GetXaxis()->GetXmax(), "MeV");
    x.setBins(p.nSamplingBins);
    x.setRange("full",x.getMin(),r.threshold+2); // for close-to-threshold tagger bins

    // load data to be fitted
    RooDataHist h_roo_data("h_roo_data","dataset",x,p.h_data);

    // build shifted mc lineshape
    RooRealVar x_shift("x_shift", "shift in IM", 0.0, -10.0, 10.0);
    RooProduct x_shift_invert("x_shift_invert","shifted IM",RooArgSet(x_shift, RooConst(-1.0)));
    RooAddition x_shifted("x_shifted","shifted IM",RooArgSet(x,x_shift_invert));
    RooDataHist h_roo_mc("h_roo_mc","MC lineshape", x, p.h_mc);
    RooHistPdf pdf_mc_lineshape("pdf_mc_lineshape","MC lineshape as PDF", x_shifted, x, h_roo_mc, p.interpOrder);

    // build detector resolution smearing

    RooRealVar  sigma("sigma","detector resolution",  2.0, 0.0, 10.0);
    RooGaussian pdf_smearing("pdf_smearing","Single Gaussian", x, RooConst(0.0), sigma);

    // build signal as convolution, note that the gaussian must be the second PDF (see documentation)
    RooFFTConvPdf pdf_signal("pdf_signal","MC_lineshape (X) gauss", x, pdf_mc_lineshape, pdf_smearing);

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
    //    RooChebychev pdf_background("chebychev","Polynomial background",x,roo_bkg_params);
    //    x.setRange("bkg_l", x.getMin(), 930);
    //    x.setRange("bkg_r", 990, 1000);
    //    pdf_background.fitTo(h_roo_data, Range("bkg_l,bkg_r"), Extended()); // using Range(..., ...) does not work here (bug in RooFit, sigh)



    RooRealVar argus_cutoff("argus_cutoff","argus pos param", r.threshold);
    RooRealVar argus_shape("argus_shape","argus shape param", -5, -25.0, 0.0);
    RooRealVar argus_p("argus_p","argus p param", 0.5);
    RooArgusBG pdf_background("pdf_background","bkg argus",x,argus_cutoff,argus_shape,argus_p);

    // build sum
    RooRealVar nsig("nsig","#signal events", 3e3, 0, 1e6);
    RooRealVar nbkg("nbkg","#background events", 3e3, 0, 1e6);
    RooAddPdf pdf_sum("pdf_sum","total sum",RooArgList(pdf_signal,pdf_background),RooArgList(nsig,nbkg));
    // do some pre-fitting to obtain better starting values, make sure function is non-zero in range
    //    x.setRange("nonzero",x.getMin(), threshold-5);
    //    pdf_sum.chi2FitTo(h_roo_data, Range("nonzero"), PrintLevel(-1)); // using Range(..., ...) does not work here (bug in RooFit, sigh)


    // do the actual maximum likelihood fit
    // use , Optimize(false), Strategy(2) for double gaussian...?!
    r.fitresult = pdf_sum.fitTo(h_roo_data, Extended(), SumW2Error(kTRUE), Range("full"), Save(), PrintLevel(-1));

    // draw output and remember pointer
    r.fitplot = x.frame();

    h_roo_data.plotOn(r.fitplot);
    r.h_data = dynamic_cast<RooHist*>(r.fitplot->findObject(0));

    // need to figure out chi2nds and stuff after plotting data and finally fitted pdf_sum
    // also the residHist must be created here (and rememebered for later use)
    pdf_sum.plotOn(r.fitplot, LineColor(kRed));
    //    pdf_sum.plotOn(frame, LineColor(kRed), VisualizeError(*fr));
    r.f_sum = dynamic_cast<RooCurve*>(r.fitplot->findObject(0));
    r.chi2ndf = r.fitplot->chiSquare(r.numParams());

    auto pdf_sum_tf = pdf_sum.asTF(RooArgList(x), RooArgList(*pdf_sum.getParameters(x)), RooArgSet(x));
    r.peakpos = pdf_sum_tf->GetMaximumX(p.signal_region.Start(), p.signal_region.Stop());
    r.residual = r.fitplot->residHist();

    pdf_sum.plotOn(r.fitplot, Components(pdf_background), LineColor(kBlue));
    r.f_bkg = dynamic_cast<RooCurve*>(r.fitplot->findObject(0));
    pdf_sum.plotOn(r.fitplot, Components(pdf_signal), LineColor(kGreen));
    r.f_sig = dynamic_cast<RooCurve*>(r.fitplot->findObject(0));
    return r;
}

template<typename T, typename Transform>
N_t calcSum(const std::vector<T>& input, Transform transform) {
    std::vector<N_t> Ns(input.size());
    std::transform(input.begin(), input.end(), Ns.begin(), transform);
    N_t N_sum; // sigma=0 means unmeasured

    APLCON::Fit_Settings_t fit_settings;
    fit_settings.ConstraintAccuracy = 1e-2;
    APLCON::Fitter<std::vector<N_t>, N_t> fitter(fit_settings);
    fitter.DoFit(Ns, N_sum, [] (const vector<N_t>& Ns, const N_t& Nsum) {
        double sum = 0.0;
        for(auto& n : Ns)
            sum += n.Value;
        return Nsum.Value - sum;
    });
    return N_sum;
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

    auto cmd_input = cmd.add<TCLAP::ValueArg<string>>("i","input","ROOT input file",true,"","rootfile");

    cmd.parse(argc, argv);
    if(cmd_verbose->isSet()) {
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());
    }
    debug = cmd_debug->isSet();

    RooMsgService::instance().setGlobalKillBelow(debug ? RooFit::WARNING : RooFit::ERROR);
    if(!debug)
        RooMsgService::instance().setSilentMode(true);

    ExpConfig::Setup::SetByName(cmd_setup->getValue());
    auto Tagger = ExpConfig::Setup::GetDetector<TaggerDetector_t>();

    const string ref_prefix   = "EtapOmegaG_plot_Ref";
    const string ref_histpath = ref_prefix+"/DiscardedEk=0/KinFitProb>0.02";
    const string ref_histname = "h_IM_2g_TaggCh";

    TH2D* ref_data;
    TH2D* ref_mc;
    TH1D* ref_mctrue_generated;

    WrapTFileInput input(cmd_input->getValue()); // keep it open
    {
        const string histpath = ref_histpath+"/h/Data/"+ref_histname;
        if(!input.GetObject(histpath, ref_data)) {
            LOG(ERROR) << "Cannot find " << histpath;
            return EXIT_FAILURE;
        }
    }
    {
        const string histpath = ref_histpath+"/h/Ref/"+ref_histname;
        if(!input.GetObject(histpath, ref_mc)) {
            LOG(ERROR) << "Cannot find " << histpath;
            return EXIT_FAILURE;
        }
    }
    {
        const string histpath = ref_prefix+"/h_mctrue_generated";
        if(!input.GetObject(histpath, ref_mctrue_generated)) {
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


    // start creating the overview (more will be added after fits)
    ant::canvas c_overview("EtapOmegaG_fit: Ref Overview");
    c_overview << drawoption("colz")
               << ref_mc << ref_data << ref_mctrue_generated
               << drawoption("");

    std::vector<fit_return_t> fit_results;

    const interval<int> taggChRange{0, 40}; // max should be 40 or 39

    for(auto taggch=taggChRange.Stop();taggch>=taggChRange.Start();taggch--) {
        LOG(INFO) << "Fitting TaggCh=" << taggch;
        fit_params_t p;
        p.TaggCh = taggch;
        p.Eg = Tagger->GetPhotonEnergy(taggch);

        // fit MC lineshape to data
        const auto taggbin = taggch+1;
        p.h_mc   = ref_mc->ProjectionX("h_mc",taggbin,taggbin);
        p.h_data = ref_data->ProjectionX("h_data",taggbin,taggbin);
        auto r = doFit_Ref(p);

        // determine efficiency corrected N_effcorr = N_fit/efficiency = N_fit * mc_generated/mc_reco;
        {
            auto N_fit = r.getNfit();
            N_t N_mcreco;
            N_mcreco.Value = p.h_mc->IntegralAndError(1, p.h_mc->GetNbinsX(), N_mcreco.Sigma, ""); // take binwidth into account?
            N_t N_mcgen(ref_mctrue_generated->GetBinContent(taggbin), ref_mctrue_generated->GetBinError(taggbin));

            APLCON::Fit_Settings_t fit_settings;
            fit_settings.ConstraintAccuracy = 1e-2;
            APLCON::Fitter<N_t, N_t, N_t, N_t> fitter(fit_settings);
            fitter.DoFit(N_fit, N_mcreco, N_mcgen, r.N_effcorr,
                         [] (const N_t& N_fit, const N_t& N_mcreco, const N_t& N_mcgen, const N_t& N_effcorr) {
                return N_effcorr.Value - N_fit.Value * N_mcgen.Value / N_mcreco.Value;
            });

            if(debug) {
                LOG(INFO) << "TaggCh=" << taggch
                          << " N_fit=" << N_fit
                          << " N_effcorr=" << r.N_effcorr;
            }
        }

        // save/plot values
        fit_results.emplace_back(r);

        if(debug) {
            LOG(INFO) << r;
        }
    }

    // plot all single fits
    ant::canvas c_plots_data("EtapOmegaG_fit: Plots Data");
    for(auto& r : fit_results)
        c_plots_data << r;
    c_plots_data << endc;

    // sum up the N_data and N_effcorr

    auto N_fit_sum = calcSum(fit_results, [] (const fit_return_t& r) {
        return r.getNfit();
    });
    LOG(INFO) << "Sum of N_fit: " << N_fit_sum;

    auto N_effcorr_sum = calcSum(fit_results, [] (const fit_return_t& r) {
        return r.N_effcorr;
    });
    LOG(INFO) << "Sum of N_fit/eff: " << N_effcorr_sum;

    N_t BR_etap_2g(2.20/100.0,0.08/100.0); // branching ratio eta'->2g is about 2.2 % (PDG)
    N_t N_etap(0,0);
    {
        APLCON::Fit_Settings_t fit_settings;
        fit_settings.ConstraintAccuracy = 1e-2;
        APLCON::Fitter<N_t, N_t, N_t> fitter(fit_settings);
        fitter.DoFit(N_effcorr_sum, BR_etap_2g, N_etap, [] (
                     const N_t& N_effcorr_sum, const N_t& BR_etap_2g, const N_t& N_etap) {
            return N_etap.Value - N_effcorr_sum.Value/BR_etap_2g.Value;
        });
    }
    LOG(INFO) << "Number of tagged eta' in EPT 2014 beamtime: " << N_etap;

    // plot summation over tagg channels, calc total chi2
    {
        const auto  makeTF1sum = [&fit_results] (RooCurve* fit_return_t::* PtrToMember) -> TF1* {
            if(fit_results.empty())
                return nullptr;
            auto& r = fit_results.front();
            auto axis = r.fitplot->GetXaxis();
            auto f = new TF1("", [&fit_results, PtrToMember] (double* x, double*) {
                double sum = 0;
                for(auto& r : fit_results) {
                    // individual thresholds of results are lower than global
                    if(x[0] < r.threshold)
                        sum += (r.*PtrToMember)->Eval(x[0]);
                }
                return sum;
            }, axis->GetXmin(), fit_results.back().threshold, 0);
            f->SetNpx(1000);
            f->SetLineColor((r.*PtrToMember)->GetLineColor());
            f->SetLineWidth((r.*PtrToMember)->GetLineWidth());
            return f;
        };

        auto& r = fit_results.front();

        auto h_proj_all_ch = ref_data->ProjectionX("proj_all_taggch",taggChRange.Start()+1,taggChRange.Stop()+1);
        h_proj_all_ch->GetYaxis()->SetTitle(r.fitplot->GetYaxis()->GetTitle());
        h_proj_all_ch->SetMarkerStyle(r.h_data->GetMarkerStyle());
        h_proj_all_ch->SetLineColor(r.h_data->GetLineColor());
        h_proj_all_ch->SetStats(0);

        auto f_sum = makeTF1sum(&fit_return_t::f_sum);

        const auto calcApproxNDF = [r, f_sum] () {
            double x_low, x_high;
            f_sum->GetRange(x_low, x_high);
            // the number of parameters for a single is used here,
            // don't know how to handle this 100% correct
            const auto numPars = r.fitresult->floatParsFinal().getSize();
            const auto binwidth = r.h_data->getNominalBinWidth();
            return (x_high-x_low)/binwidth - numPars;
        };

        const auto chi2 = h_proj_all_ch->Chisquare(f_sum,"R");
        const auto ndf = calcApproxNDF();
        auto lbl = new TPaveText();
        lbl->SetX1NDC(0.58);
        lbl->SetX2NDC(0.88);
        lbl->SetY1NDC(0.6);
        lbl->SetY2NDC(0.85);
        lbl->SetBorderSize(0);
        lbl->SetFillColor(kWhite);
        lbl->AddText(static_cast<string>(std_ext::formatter() << setprecision(0) << fixed << "N = " << N_fit_sum.Value << " #pm " << N_fit_sum.Sigma).c_str());
        lbl->AddText(static_cast<string>(std_ext::formatter() << setprecision(0) << fixed << "N/#varepsilon = " << N_effcorr_sum.Value << " #pm " << N_effcorr_sum.Sigma).c_str());
        lbl->AddText(static_cast<string>(std_ext::formatter() << setprecision(0) << fixed << "#chi^{2}_{red} = " << chi2 << "/" << ndf << " #approx " << setprecision(2) << chi2/ndf).c_str());

        c_overview << drawoption("E1") << h_proj_all_ch
                   << samepad << makeTF1sum(&fit_return_t::f_sig)
                   << samepad << makeTF1sum(&fit_return_t::f_bkg)
                   << samepad << f_sum
                   << samepad << lbl;

    }

    // number of events and effcorr events as multigraph
    {
        auto g_N_fit = new TGraphErrors();
        g_N_fit->SetTitle("N");
        g_N_fit->SetFillColor(kWhite);
        g_N_fit->SetLineColor(kRed);
        g_N_fit->SetLineWidth(3);

        auto g_N_effcorr = new TGraphErrors();
        g_N_effcorr->SetTitle("N/#varepsilon");
        g_N_effcorr->SetFillColor(kWhite);
        g_N_effcorr->SetLineColor(kBlue);
        g_N_effcorr->SetLineWidth(3);

        for(const fit_return_t& r : fit_results) {
            auto n = g_N_fit->GetN();
            auto N_fit = r.getNfit();
            g_N_fit->SetPoint(n, r.p.Eg, N_fit.Value);
            g_N_fit->SetPointError(n, Tagger-> GetPhotonEnergyWidth(r.p.TaggCh)/2, N_fit.Sigma);
            g_N_effcorr->SetPoint(n, r.p.Eg, r.N_effcorr.Value);
            g_N_effcorr->SetPointError(n, Tagger-> GetPhotonEnergyWidth(r.p.TaggCh)/2, r.N_effcorr.Sigma);
        }

        auto multigraph = new TMultiGraph("g_N_fit_effcorr","Ref: Number of events");
        multigraph->Add(g_N_fit);
        multigraph->Add(g_N_effcorr);

        c_overview << drawoption("AP") << padoption::Legend << multigraph << endc;
        // change axis properties after drawing (before the axes don't exist in ROOT...sigh)
        multigraph->GetXaxis()->SetTitle("E_{#gamma} / MeV");
        multigraph->GetYaxis()->SetTitle("Events");
        multigraph->SetMinimum(0);
        multigraph->SetMaximum(4400);
        // necesary to immediately show changes to multigraph after drawing in canvas
        gPad->Modified();
        gPad->Update();
    }

    // run TRint
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

    return EXIT_SUCCESS;
}
