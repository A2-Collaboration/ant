#include "analysis/physics/etaprime/etaprime_dalitz.h"

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cassert>
#include <fstream>
#include <filesystem>

#include "APLCON.hpp"

#include "base/Logger.h"

#include "tclap/CmdLine.h"
#include "tclap/ValuesConstraintExtra.h"
#include "base/interval.h"
#include "base/WrapTFile.h"
#include "base/std_ext/string.h"
#include "base/std_ext/system.h"
#include "base/std_ext/math.h"
#include "base/ParticleType.h"

#include "analysis/plot/RootDraw.h"
#include "analysis/plot/HistogramFactory.h"
#include "analysis/utils/ParticleTools.h"
#include "expconfig/ExpConfig.h"
#include "base/Detector_t.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"

#include "TSystem.h"
#include "TRint.h"
#include "TROOT.h"
#include "TStyle.h"

#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooGaussian.h"
#include "RooArgusBG.h"
#include "RooCBShape.h"
#include "RooNovosibirsk.h"
#include "root-addons/roofit_extensions/RooGaussExp.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooDataHist.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "RooFFTConvPdf.h"
#include "RooChi2Var.h"
#include "RooMinuit.h"
#include "RooFitResult.h"

#include "detail/tools.h"


namespace fs = std::filesystem;

using namespace ant;
using namespace std;
using namespace RooFit;

using q2_params_t = ant::analysis::physics::EtapDalitzTools::q2_params_t;


static volatile sig_atomic_t interrupt = false;


class Fitter {

public:

    explicit Fitter(const APLCON::Fit_Settings_t& settings = {}) :
        aplcon(settings)
    {}

    // Value, Sigma, Pull
    struct N_etap_t {
        explicit N_etap_t(double v = 0, double s = 0) :
            Value(v), Sigma(s) {}

        double Value;
        double Sigma = std_ext::NaN;
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

        void SetValueSigma(double value, double sigma) {
            Value = value;
            Sigma = sigma;
        }

        operator const double&() const noexcept {
            return Value;
        }

        friend ostream& operator<<(ostream& s, const N_etap_t& n) {
            return s << n.Value << " +/- " << n.Sigma;
        }
    };

    APLCON::Fitter<vector<N_etap_t>, N_etap_t> aplcon;  // template parameters: list of individual fits, combined fit result

    APLCON::Result_t DoFit(vector<N_etap_t>& N, N_etap_t& sum)
    {
        const auto& r = aplcon.DoFit(N, sum, [] (const vector<N_etap_t>& N, const N_etap_t& N_sum) {
            double sum = 0.;
            for (const auto& n : N)
                sum += n.Value;
            return N_sum.Value - sum;
        });
        return r;
    }
};

// C++11 compliant thread-safe singleton
struct settings_t final
{
    static settings_t& get();

    bool debug;
    int rebin;
    fs::path out_dir;

private:
    settings_t() = default;
    ~settings_t() = default;

    // delete copy and move constructors
    settings_t(const settings_t&) = delete;
    settings_t& operator=(const settings_t&) = delete;
    settings_t(settings_t&&) = delete;
    settings_t& operator=(settings_t&&) = delete;
};

settings_t& settings_t::get()
{
    // create a magic static
    static settings_t instance;
    return instance;
}


// convenience method to print vectors
template<typename T>
ostream& operator<< (ostream& out, const vector<T>& v)
{
    out << "{";
    auto it = begin(v);
    for (const auto& i : v) {
        out << i;
        if (++it != v.end())
            out << ", ";
    }
    out << "}";
    return out;
}

template <typename T>
vector<T> convert_piecewise_interval(const PiecewiseInterval<T>& interval, const bool make_unique = true)
{
    vector<T> v;
    for (const auto& range : interval)
        for (auto i = range.Start(); i <= range.Stop(); i++)
            v.emplace_back(i);

    if (!make_unique)
        return v;

    // sort the vector and remove duplicate entries
    sort(v.begin(), v.end());
    auto it = unique(v.begin(), v.end());
    auto dist = unsigned(distance(v.begin(), it));
    auto duplicates = v.size() - dist;
    if (duplicates) {
        LOG(WARNING) << "The provided intervals contain " << duplicates << " duplicate entries, they will be removed";
        v.resize(dist);
    }

    return v;
}

template <typename T, typename U>
vector<T> convert_piecewise_interval_type(const PiecewiseInterval<U>& interval, const bool make_unique = true)
{
    const auto v = convert_piecewise_interval(interval, make_unique);
    vector<T> vt;
    transform(v.begin(), v.end(), back_inserter(vt), [] (U u) { return static_cast<T>(u); });
    return vt;
}


string concat_string(const vector<string>& strings, const string& delimiter = ", ")
{
    if (strings.empty())
        return "";

    return accumulate(next(strings.begin()), strings.end(), strings.front(),
            [&delimiter] (string& concat_str, const string& str) {
                return concat_str + delimiter + str;
            });
}

string cuts_path(const vector<string>& cuts, const char* delimiter = "/")
{
    if (cuts.empty())
        return "";

    stringstream s;
    copy(cuts.begin(), cuts.end(), ostream_iterator<string>(s, delimiter));
    return s.str();
}

string get_path(const string& cut_string, const string& tree)
{
    return tree + "/" + cut_string;
}


struct q2_bin_cut_t {
    string q2_bin;
    vector<string> cuts;
    string cut_string;

    void create_cut_string() {
        cut_string = cuts_path(cuts, "/");
    }

    q2_bin_cut_t(const string& bin, const string& cuts) : q2_bin(bin), cut_string(cuts) {}

    q2_bin_cut_t(const string& bin, const vector<string>& _cuts) : q2_bin(bin), cuts(_cuts) {
        create_cut_string();
    }
};

struct fit_result_t {
    int taggCh = -1;
    double chi2ndf = std_ext::NaN;
    double n_etap, n_error, eff_corr;

    RooCurve* signal = nullptr;
    RooCurve* bkg = nullptr;
};


void test_path_building()
{
    const vector<string> cuts = {
        "selection",
        "KinFitProb > 0.1",
        "nothing",
        "thight cluster size"};
    const string tree = "EtapDalitz_plot_Sig";
    cout << "Test building tree path from cuts vector:" << endl
        << get_path(concat_string(cuts, "/"), tree) << endl;
    cout << "Test cuts path using copy and stringstream:" << endl
        << cuts_path(cuts) << endl;
}

void traverseCuts(TDirectory* dir, vector<vector<string>>& cuts) {
    auto keys = dir->GetListOfKeys();
    if (!keys)
        return;

    vector<string> dirnames;
    bool h_found = false;
    TIter nextk(keys);
    TKey* key;
    TKey* nextdir = nullptr;
    while ((key = static_cast<TKey*>(nextk())))
    {
        auto classPtr = TClass::GetClass(key->GetClassName());
        if (classPtr->InheritsFrom(TDirectory::Class())) {
            const string dirname(key->GetName());
            if (dirname == "h")
                h_found = true;
            else {
                nextdir = key;
                dirnames.emplace_back(dirname);
            }
        }
    }

    if (h_found && !dirnames.empty()) {
        cuts.emplace_back(dirnames);
        if (nextdir) {
            traverseCuts(dynamic_cast<TDirectory*>(nextdir->ReadObj()), cuts);
        }
    }
}

vector<vector<string>> extractCuts(const string& prefix, const WrapTFileInput& input) {
    TDirectory* prefixDir = nullptr;
    if (!input.GetObject(prefix, prefixDir))
        throw runtime_error("Cannot find prefix dir " + prefix);
    vector<vector<string>> cuts;
    traverseCuts(prefixDir, cuts);
    return cuts;
}

void print_extracted_cuts(const string& file)
{
    WrapTFileInput input(file);
    const string prefix = "EtapDalitz_plot_Sig";
    auto cuts = extractCuts(prefix, input);
    cout << "Extracted Cuts:" << endl;
    size_t cut_level = 0;
    for (const auto& vec : cuts) {
        cout << "  cut level " << ++cut_level << endl;
        for (const auto& cut : vec)
            cout << "    " << cut << endl;
    }
}



template <typename T>
struct draw_TGraph_t : ant::root_drawable_traits {
    T* graph;
    string xlabel;
    string ylabel;
    interval<double> yrange;

    explicit draw_TGraph_t(T* g, const string& xlabel_, const string& ylabel_ = "",
                           const interval<double>& yrange_ = {0,-1}) :
        graph(g), xlabel(xlabel_), ylabel(ylabel_), yrange(yrange_)
    {}

    void Draw(const string& opt) const override
    {
        graph->Draw(opt.c_str());
        graph->GetXaxis()->SetTitle(xlabel.c_str());
        graph->GetYaxis()->SetTitle(ylabel.c_str());
        if (yrange.IsSane()) {
            graph->SetMinimum(yrange.Start());
            graph->SetMaximum(yrange.Stop());
        }
        // necessary to immediately show changes to multigraph after drawing in canvas
        gPad->Modified();
        gPad->Update();
    }
};

template <typename T, typename... Args>
draw_TGraph_t<T> draw_TGraph(T* g, Args&&... args) {
    return draw_TGraph_t<T>(g, std::forward<Args>(args)...);
}

// helper method to save a TCanvas, TPad, TVirtualPad, ... to a file if the directory exists
template <typename T>
void save_pad(const T* const p, const fs::path& output_dir, const string& filename)
{
    if (output_dir.empty())
        return;

    const auto path = output_dir / fs::path(filename);
    p->Print(path.c_str());
}

void reference_fit(const WrapTFileInput& input, const string& cuts, const vector<int>& EPTrange,
                   const WrapTFileInput& mc, pair<double, double>& result)
{
    settings_t& settings = settings_t::get();

    TH2D* ref_data;
    TH2D* ref_mc;
    TH1* h_data;
    TH1* h_mc;
    TH2D* trueIM_EPT = nullptr;
    TH1* h_true = nullptr;

    // check if MC file provided, if yes get true MC histogram for EPT vs IM
    if (mc.NumberOfFiles()) {
        if (!mc.GetObject("Etap2gMC/h_taggCh_vs_trueIM", trueIM_EPT))
            throw runtime_error("Couldn't find true MC histogram in file " + mc.FileNames());
    } else
        LOG(WARNING) << "No MC input provided, some default values will be used for efficiency corrections";

    TCanvas* c = new TCanvas("c", "", 10,10, 800,800);

    string hist = "EtapDalitz_plot_Ref/" + cuts +  "/h/Data/taggChannel_vs_etapIM_kinfitted";
    if (!input.GetObject(hist, ref_data))
        throw runtime_error("Couldn't find " + hist + " in file " + input.FileNames());

    hist = "EtapDalitz_plot_Ref/" + cuts +  "/h/Reference/taggChannel_vs_etapIM_kinfitted";
    if (!input.GetObject(hist, ref_mc))
        throw runtime_error("Couldn't find " + hist + " in file " + input.FileNames());

    auto EPT = ExpConfig::Setup::GetDetector<TaggerDetector_t>();

    const auto maxIM = [] (const double Eg) {
        const auto mp = ParticleTypeDatabase::Proton.Mass();
        return sqrt(mp*mp + 2*mp*Eg) - mp;
    };

    constexpr IntervalD fit_range = {840, 1020};

    double total_number_etap = 0.;
    double total_n_err = 0.;

    vector<fit_result_t> results;

    auto g_n = new TGraphErrors();
    g_n->SetTitle("Number #eta'");
    g_n->SetLineColor(kRed);
    g_n->SetLineWidth(2);
    g_n->SetMarkerSize(0);

    canvas c_N("Number eta' based on Reference");

    // loop over the provided EPT channels
    for (const auto taggCh : EPTrange) {
        if (interrupt)
            break;

        fit_result_t res;
        res.taggCh = taggCh;

        const double taggE = EPT->GetPhotonEnergy(unsigned(taggCh));
        const int taggBin = taggCh+1;
        LOG(INFO) << "Fitting EPT channel " << taggCh << " (E_gamma = " << taggE << " MeV)";

        h_data = ref_data->ProjectionX("h_data", taggBin, taggBin);
        if (taggCh == 40)  // close to threshold, decrease histogram IM range
            h_data->GetXaxis()->SetRangeUser(900,1100);
        h_mc = ref_mc->ProjectionX("h_mc", taggBin, taggBin);
        if (trueIM_EPT)
            h_true = trueIM_EPT->ProjectionX("h_true", taggBin, taggBin);

        const double cutoff = maxIM(taggE);
        VLOG(1) << "EPT E = " << taggE << ", calculated cutoff value: " << cutoff;

        // clear and prepare canvas
        c->Clear();
        c->SetTitle(Form("Fit: %s", h_data->GetTitle()));
        c->cd();
        c->Divide(1,2);
        c->cd(1);

        // define observable and ranges
        RooRealVar var_IM("IM","IM", fit_range.Start(), fit_range.Stop(), "MeV");
        var_IM.setBins(1000);
        var_IM.setRange("full", fit_range.Start(), fit_range.Stop());  // define "full" range used for fitting

        // load data to be fitted
        if (settings.rebin)
            h_data->Rebin(settings.rebin);
        RooDataHist h_roo_data("h_roo_data","dataset",var_IM,h_data);

        /* MC lineshape fit
        // build shifted mc lineshape
        const double offset = h_data->GetBinCenter(h_data->GetMaximumBin()) - ParticleTypeDatabase::EtaPrime.Mass();
        RooRealVar var_IM_shift("var_IM_shift", "shift in IM", offset, -20., 20.);  // use current offset as starting value (just using 0 would work equally fine)
        RooProduct var_IM_shift_invert("var_IM_shift_invert","shifted IM",RooArgSet(var_IM_shift, RooConst(-1.)));
        RooAddition var_IM_shifted("var_IM_shifted","shifted IM",RooArgSet(var_IM,var_IM_shift_invert));
        RooDataHist h_roo_mc("h_roo_mc","MC lineshape", var_IM, h_mc);
        RooHistPdf pdf_mc_lineshape("pdf_mc_lineshape","MC lineshape as PDF", var_IM_shifted, var_IM, h_roo_mc, 2);  // 2nd order interpolation (or 4th?)

        // build gaussian
        RooRealVar  var_gauss_sigma("gauss_sigma","detector resolution", 5., .01, 20.);
        RooGaussian pdf_gaussian("pdf_gaussian","Gaussian smearing", var_IM, RooConst(0.), var_gauss_sigma);

        // build signal as convolution, note that the gaussian must be the second PDF (see documentation)
        RooFFTConvPdf pdf_signal("pdf_signal","MC_lineshape (X) gauss",var_IM, pdf_mc_lineshape, pdf_gaussian) ;
        */

        // --- Build CB Function PDF ---
        RooRealVar cb_x0("cb_x0", "expectation value", 960, 950, 975);
        RooRealVar cb_sigma("cb_sigma", "standard deviation", 3, .01, 20);
        RooRealVar cb_alpha("cb_alpha", "transition gauss to power function", 1.3, .5, 2.);
        RooRealVar cb_n("cb_n", "parameter power function", 1, .1, 10);
        //RooCBShape pdf_signal("signal", "CB Function", var_IM, cb_x0, cb_sigma, cb_alpha, cb_n);
        //RooNovosibirsk pdf_signal("signal", "Novosibirsk function", var_IM, cb_x0, cb_sigma, cb_n);
        RooGaussExp pdf_signal("signal", "Simple CB function", var_IM, cb_x0, cb_sigma, cb_n);

        // build background with ARGUS function
        RooRealVar argus_cutoff("argus_cutoff","argus pos param", cutoff);  // upper threshold, calculated for beam energy
        RooRealVar argus_shape("argus_chi","argus shape param #chi", -5, -25., 5.);
        //RooRealVar argus_p("argus_p","argus p param", 0.5, 0, 1);
        RooRealVar argus_p("argus_p","argus p param", .5);
        RooArgusBG pdf_background("pdf_background","bkg argus",var_IM,argus_cutoff,argus_shape,argus_p);

        const double n_total = h_data->Integral();
        // build sum
        RooRealVar nsig("N_sig","#signal events", n_total/2, 0., 2*n_total);
        RooRealVar nbkg("N_bkg","#background events", n_total/2, 0., 2*n_total);
        RooAddPdf pdf_sum("pdf_sum","total sum",RooArgList(pdf_signal,pdf_background),RooArgList(nsig,nbkg));

        RooFitResult* fit = pdf_sum.fitTo(h_roo_data, Extended(), SumW2Error(kTRUE), Range("full"), Save(), PrintLevel(-1));
        fit->Print();
        LOG(INFO) << "Covariance Matrix:";
        fit->covarianceMatrix().Print();
        RooPlot* frame = var_IM.frame();
        h_roo_data.plotOn(frame);
        frame->GetXaxis()->SetLabelSize(.05f);
        frame->GetXaxis()->SetTitleSize(.05f);
        frame->GetYaxis()->SetLabelSize(.05f);
        frame->GetYaxis()->SetTitleSize(.05f);
        frame->GetXaxis()->SetRangeUser(fit_range.Start(), fit_range.Stop());
        frame->SetTitle("Reference");

        auto p = new TPaveText();
        p->SetFillColor(0);
        p->SetFillStyle(0);
        p->SetX1NDC(0.14);
        p->SetX2NDC(0.39);
        p->SetY1NDC(0.38);
        p->SetY2NDC(0.86);
        p->SetTextSize(.04f);

        // define lambda to insert lines in stat box
        const auto addLine = [] (TPaveText& p, const RooRealVar& v, const string& name = "") {
            p.InsertText(Form("%s = %.2f #pm %.2f", name.empty() ? v.GetName() : name.c_str(), v.getValV(), v.getError()));
        };

        //pdf_background.plotOn(frame);
        pdf_sum.plotOn(frame, Components(pdf_background), Name("bkg"), LineColor(kAzure-3), PrintEvalErrors(-1));
        pdf_sum.plotOn(frame, Components(pdf_signal), Name("signal"), LineColor(kGreen+1));
        pdf_sum.plotOn(frame, Name("sum"), LineColor(kRed+1), PrintEvalErrors(-1));
        frame->Draw();
        pdf_sum.paramOn(frame);
        // number of item index seems to reflect the order plotOn is called on a RooPlot:
        // index 0 is histogram, 1 is background, 2 is signal, 3 is sum (called ploton in that order), 4 is paramBox
        res.signal = frame->getCurve("signal");
        res.bkg = frame->getCurve("bkg");

        RooHist* hresid = frame->residHist();
        hresid->SetTitle("Residuals");
        hresid->GetXaxis()->SetRangeUser(fit_range.Start(), fit_range.Stop());
        hresid->GetXaxis()->SetTitle("m(#gamma#gamma) [MeV]");
        hresid->GetXaxis()->SetLabelSize(.05f);
        hresid->GetXaxis()->SetTitleSize(.05f);
        hresid->GetXaxis()->SetTickLength(.08f);
        hresid->GetYaxis()->SetLabelSize(.05f);

        const double chi2ndf = frame->chiSquare(fit->floatParsFinal().getSize());
        res.chi2ndf = chi2ndf;

        p->InsertText(Form("#chi^{2}/dof = %.2f", chi2ndf));
//        addLine(*p, var_IM_shift,    "#Delta IM");
//        addLine(*p, var_gauss_sigma, "#sigma");
        addLine(*p, cb_x0,    "cb_mean");
        addLine(*p, cb_sigma, "cb_sigma");
        addLine(*p, cb_alpha, "cb_alpha");
        addLine(*p, cb_n,     "cb_n");
        addLine(*p, argus_cutoff,    "c");
        addLine(*p, argus_shape,     "#chi");
        addLine(*p, argus_p,         "p");
        addLine(*p, nsig,            "n_{sig}");
        addLine(*p, nbkg,            "n_{bkg}");
        p->Draw();

        c->cd(2);
        hresid->Draw();

        //TODO: maybe provide default efficiency corrections if no true MC histogram provided
        constexpr double BR2g = .022;
        const double eff_corr = h_true ? h_mc->GetEntries()/h_true->GetEntries() : 35034360./1e8;
        const double n_tot_corr = nsig.getValV()/eff_corr/BR2g;
        const double n_error = nsig.getError()/eff_corr/BR2g;
        total_number_etap += n_tot_corr;
        total_n_err += n_error*n_error;
        res.n_etap = n_tot_corr;
        res.n_error = n_error;
        res.eff_corr = eff_corr;
        LOG(INFO) << "Number of efficiency corrected eta' for EPT channel "
                  << taggCh << ": " << n_tot_corr << " +/- " << n_error;

        c->Modified();
        c->Update();

        save_pad(c, settings.out_dir, "ref_fit_channel" + to_string(taggCh) + ".pdf");

        // add the number of eta' for the current EPT channel to the corresponding graph
        // clear the canvas to update the plotted graph for each fit within the loop
        {
            const int n = g_n->GetN();
            g_n->SetPoint(n, taggE, n_tot_corr);
            g_n->SetPointError(n, EPT->GetPhotonEnergyWidth(unsigned(taggCh))/2., n_error);
            c_N.clear();
            c_N << drawoption("AP") << draw_TGraph(g_n, "E_{#gamma} [MeV]", "##eta' / EPT Ch.", interval<double>{0,21e4}) << endc;
        }

        results.emplace_back(move(res));
    }

    result = {total_number_etap, sqrt(total_n_err)};
    LOG(INFO) << "Total number of eta': " << result.first << " +/- " << result.second;

    c_N.cd();
    save_pad(gPad, settings.out_dir, "ref_Netap_vs_Eg.pdf");
    // write canvas with number of eta' per energy bin to the output file (needed because it's a ant::canvas)
    gPad->Write("N_etap_vs_chE");


    hist = "EtapDalitz_plot_Ref/" + cuts +  "/h/Data/etapIM_kinfitted";
    if (!input.GetObject(hist, h_data))
        throw runtime_error("Couldn't find " + hist + " in file " + input.FileNames());

    // sum up all signal and background fits
    RooCurve sigSum = *results.front().signal;
    RooCurve bgSum = *results.front().bkg;
    const auto n_res = results.size();
    for (unsigned i = 1; i < n_res-1; i++) {
        sigSum = RooCurve("", "", sigSum, *results.at(i).signal);
        bgSum = RooCurve("", "", bgSum, *results.at(i).bkg);
    }
    RooCurve* signalSum = new RooCurve("signalSum", "Sum of Signal Curves", sigSum, *results.back().signal);
    RooCurve* bkgSum = new RooCurve("bkgSum", "Sum of Background Curves", bgSum, *results.back().bkg);
    RooCurve* sum = new RooCurve("sum", "Sum of all Background and Signal Curves", *signalSum, *bkgSum);

    RooRealVar var_IM("IM","IM", fit_range.Start(), fit_range.Stop(), "MeV");
    var_IM.setBins(1000);
    var_IM.setRange("full", fit_range.Start(), fit_range.Stop());

    // load data
    if (settings.rebin)
        h_data->Rebin(settings.rebin);
    RooDataHist h_roo_data("h_roo_data","dataset",var_IM,h_data);

    RooPlot* frame = var_IM.frame();
    h_roo_data.plotOn(frame);
    sum->SetLineColor(kRed+1);
    bkgSum->SetLineStyle(kDashed);
    frame->addObject(bkgSum);
    frame->addObject(sum);

    frame->GetXaxis()->SetLabelSize(.05f);
    frame->GetXaxis()->SetTitleSize(.05f);
    frame->GetYaxis()->SetLabelSize(.05f);
    frame->GetYaxis()->SetTitleSize(.05f);
    frame->GetYaxis()->SetTitleOffset(1.25f);
    frame->SetTitle("All EPT Channel Fits Combined");

    TCanvas *cSum = new TCanvas("cSum", "Sum test", 850, 10, 1000, 800);
    cSum->SetLeftMargin(0.13f);
    frame->Draw();
    cSum->Update();

    save_pad(cSum, settings.out_dir, "ref_fit_sum.pdf");


    if (interrupt)
        return;

    using N_etap = Fitter::N_etap_t;
    APLCON::Fit_Settings_t fit_settings;
    fit_settings.ConstraintAccuracy = 1e-4;
    Fitter fit(fit_settings);
    vector<N_etap> N(results.size());
    transform(results.begin(), results.end(), N.begin(), [] (const fit_result_t& r) {
        return N_etap(r.n_etap, r.n_error);
    });
    N_etap N_fitted;
    fit.DoFit(N, N_fitted);
    const auto& aplcon_res = fit.DoFit(N, N_fitted);

//TODO: status seems not to be success, but result looks reasonable --> understand what result is and why
//    if (aplcon_res.Status != APLCON::Result_Status_t::Success)
//        LOG(FATAL) << "Fit didn't work, combining individual fits failed";

    LOG(INFO) << "Fitted N eta' via APLCON: " << N_fitted;
    VLOG(1) << "Used Iterations: " << aplcon_res.NIterations << "; Fit Probability: " << aplcon_res.Probability;
}

vector<IntervalD> get_q2_ranges()
{
    vector<IntervalD> q2_ranges;

    double bin_start = q2_params_t::min_value;
    auto it = q2_params_t::bin_widths.begin();

    while (bin_start < q2_params_t::max_value) {
        // sanity check to make sure enough bin widths are provided to cover the whole region
        if (it == q2_params_t::bin_widths.end()) {
            LOG(ERROR) << "Not enough bins provided, max q^2 value of " << q2_params_t::max_value << " not reached";
            LOG(INFO) << "Given bin widths only cover the region up until " << bin_start;
            throw runtime_error("Not enough bins provided");
        }

        double q2 = bin_start + *it++;
        q2_ranges.emplace_back(bin_start, q2);

        bin_start = q2;
    }

    return q2_ranges;
}

vector<string> build_q2_histnames()
{
    vector<IntervalD> q2_ranges = get_q2_ranges();
    vector<string> hist_names;

    for (auto range : q2_ranges) {
        stringstream name;
        name << "imee_" << range.Start() << "_" << range.Stop();
        hist_names.emplace_back(name.str());
    }

    return hist_names;
}

// some static declarations to control the signal fits, not ideal this way
// TODO: rewrite this in a proper way without globally defined static objects
static constexpr IntervalD exclude_range = {950, 970};
static constexpr IntervalD default_fit_range = {840, 1020};
static TF1* signal_bg = nullptr;

double bkg_fun(double* x, double* par)
{
    if (!signal_bg)
        signal_bg = new TF1("bkg_function", "pol3", default_fit_range.Start(), default_fit_range.Stop());

    if (x[0] > exclude_range.Start() && x[0] < exclude_range.Stop()) {
        TF1::RejectPoint();
        return 0;
   }

    signal_bg->SetParameters(par);

    return signal_bg->Eval(x[0]);
}

double positive_pol(double* x, double* par)
{
    signal_bg->SetParameters(par);
    return signal_bg->Eval(x[0]) < 0 ? 0 : signal_bg->Eval(x[0]);
}

double normed_gauss(double* x, double* par)
{
    static const double sqrt2pi = sqrt(2*M_PI);

    const double exponent = (x[0]-par[1])/par[2];
    return par[0]/par[2]/sqrt2pi*exp(-.5*exponent*exponent);
}

double normed_gauss_binwidth(double* x, double* par)
{
    static const double sqrt2pi = sqrt(2*M_PI);

    const double exponent = (x[0]-par[1])/par[2];
    return par[0]*par[3]/par[2]/sqrt2pi*exp(-.5*exponent*exponent);
}

double total(double* x, double* par)
{
    return positive_pol(x, par) + normed_gauss_binwidth(x, &par[4]);
}

TH1* subtract_bkg(const TH1* const h, TF1* const f)
{
    TH1* hist = dynamic_cast<TH1D*>(h->Clone("bkg_subtracted"));
    // IMPORTANT!!! Calling this prevents the histogram from being empty,
    // although documentation says otherwise, it's still ROOT we're talking about... -.-
    hist->SetDirectory(nullptr);

    bool res = hist->Add(f, -1.);
    if (!res) {
        LOG(ERROR) << "Something went wrong while subtracting the background...";
        return nullptr;
    }
    return hist;
}

void signal_fit(const WrapTFileInput& input, const vector<vector<string>>& cuts, const vector<unsigned>& imee_bins,
                const WrapTFileInput& mc, vector<pair<double, double>>& signal_fits)
{
    settings_t& settings = settings_t::get();

    auto hist_names = build_q2_histnames();
    LOG(INFO) << "size: " << hist_names.size() << "; contents:  " << hist_names;

    const auto debug = el::Loggers::verboseLevel();
    if (debug)
        cout << "Constructed q2 histogram names: " << hist_names << endl;

    if (cuts.empty())
        LOG(INFO) << "Use default cuts";

    const vector<string> cuts_q2bins = {
        "selection/KinFitProb > 0.05/cluster size/free vz cut 6/PID time cut/treefit vz cut tighter",  // 20-70
        "selection/KinFitProb > 0.05/cluster size/free vz cut 6/PID time cut/treefit vz cut tighter",  // 70-110
        "selection/KinFitProb > 0.05/cluster size/free vz cut 6/PID time cut/treefit vz cut tighter",  // 110-150
        "selection/KinFitProb > 0.05/cluster size/free vz cut 4/PID time cut/treefit vz cut tighter",  // 150-200
        "selection/KinFitProb > 0.05/cluster size/free vz cut 4/PID time cut/treefit vz cut tighter",  // 200-250
        "selection/KinFitProb > 0.05/cluster size/free vz cut 4/PID time cut/treefit vz cut tighter",  // 250-300
        "selection/KinFitProb > 0.1/cluster size/free vz cut 6/PID time cut/treefit vz cut tighter",   // 300-360
        "selection/KinFitProb > 0.1/cluster size/free vz cut 4/PID time cut/treefit vz cut tighter",   // 360-430
        "selection/KinFitProb > 0.15/tight cluster size/free vz cut 4/PID time cut/treefit vz cut tighter",  // 430-510
        "selection/KinFitProb > 0.2/tight cluster size/free vz cut 4/PID lepton cut tighter",  // 510-620
        "selection/KinFitProb > 0.2/tight cluster size/free vz cut 4/PID lepton cut tighter",  // 620-700
        "selection/KinFitProb > 0.2/tight cluster size/free vz cut 4/PID time cut/treefit vz cut tighter",  // 700-750
        "selection/KinFitProb > 0.2/tight cluster size/free vz cut 4/PID time cut/bigger lateral moment",  // 750-800
        "selection/KinFitProb > 0.2/tight cluster size/free vz cut 4/PID time cut/bigger lateral moment"  // 800-860
    };

    if (hist_names.size() < cuts_q2bins.size())
        throw runtime_error("More cut strings constructed than q2 bins");
    if (imee_bins.back() >= cuts_q2bins.size())
        throw runtime_error("Largest q2 bin index provided is too big (" + to_string(imee_bins.back())
                            + " given, but size of used cut strings vector is just " + to_string(cuts_q2bins.size()) + ")");

    vector<q2_bin_cut_t> q2_bins_cuts;
    for (auto bin = hist_names.cbegin(), cut = cuts_q2bins.cbegin();
         bin != hist_names.cend() && cut != cuts_q2bins.cend(); ++bin, ++cut)
        q2_bins_cuts.emplace_back(*bin, *cut);

    for (auto i : q2_bins_cuts)
        cout << "hist " << i.q2_bin << " will use the cut " << i.cut_string << endl;

    if (debug)
        for (const auto& i : imee_bins)
            cout << "use bin " << i << " (" << q2_bins_cuts.at(i).q2_bin << ")" << endl;

    const vector<IntervalD> fit_range = {
        {840, 1020},
        {840, 1020},
        {870, 1010},  //  2: 110-150
        {910, 1020},  //  3: 150-200
        {850, 1020},  //  4: 200-250
        {840, 1020},
        {840, 1020},
        {890, 1020},  //  7: 360-430
        {910, 1020},  //  8: 430-510
        {910, 1020},
        {920, 1020},  // 10: 620-700
        {910, 1020},
        {920, 1010},  // 12: 750-800
        {840, 1020}
    };

    const vector<IntervalD> signal_fit = {
        {920, 1000},
        {920, 1000},
        {940, 1000},  //  2: 110-150
        {940, 1000},  //  3: 150-200
        {940, 1000},  //  4: 200-250
        {940, 1000},
        {920, 1000},
        {930, 1000},  //  7: 360-430
        {920, 1000},
        {940, 1000},  //  9: 510-620
        {920, 1000},
        {930, 1000},  // 11: 700-750
        {940, 1000},
        {920, 1000}
    };

    vector<pair<double, double>> imee_fits;


    const string tree_name = "EtapDalitz_plot_Sig";

    constexpr IntervalD combined_plot_xRange = {820, 1050};

    TCanvas *c1 = new TCanvas("c1", "Background Fit", 10, 10, 1000, 800);
    TCanvas *c2 = new TCanvas("c2", "Background Subtraction", 1010, 10, 1000, 800);

    const char* default_fit_options = debug ? "R0" : "R0Q";  // R use specified function range, 0 do not plot fit result, Q quiet mode

    for (auto bin : imee_bins) {

    //constexpr size_t bin = 13;
    // used fit ranges
    const double fit_min = fit_range.at(bin).Start();
    const double fit_max = fit_range.at(bin).Stop();
    const double signal_min = signal_fit.at(bin).Start();
    const double signal_max = signal_fit.at(bin).Stop();

    if (signal_bg)
        delete signal_bg;
    signal_bg = new TF1("bkg_function", "pol3", fit_min, fit_max);


    TH1* h_data;
    TH1* h_signal;

    auto cut = q2_bins_cuts.at(bin).cut_string;
    auto q2_hist = q2_bins_cuts.at(bin).q2_bin;
    auto hist = concat_string({tree_name, cut, "h/Data", q2_hist}, "/");
    if (!input.GetObject(hist, h_data))
        throw runtime_error("Couldn't find " + hist + " in file " + input.FileNames());

    hist = concat_string({tree_name, cut, "h/Signal", q2_hist}, "/");
    if (!input.GetObject(hist, h_signal))
        throw runtime_error("Couldn't find " + hist + " in file " + input.FileNames());

    // clear and prepare canvases
    c1->Clear();
    c2->Clear();
    c1->SetTitle(Form("Background Fit: %s", q2_hist.c_str()));
    c2->SetTitle(Form("Background Subtraction: %s", q2_hist.c_str()));
    c1->cd();

    h_data->Rebin(2);
    h_data->SetDirectory(nullptr);
    h_signal->SetDirectory(nullptr);

    TF1* bkg = new TF1("bkg", bkg_fun, fit_min, fit_max, 4);


    h_data->Fit(bkg, default_fit_options);
    h_data->Fit(bkg, Form("%sL", default_fit_options));  // L likelihood fit (fails less often if regular fit is run before); WL: weighted likelihood
    TF1* bg = new TF1("bg", positive_pol, fit_min, fit_max, 4);
    bg->SetParameters(bkg->GetParameters());

    h_data->Draw();
    bg->SetLineColor(kBlue);
    bg->Draw("SAME");
    bkg->Draw("SAME");
    c1->Update();
    save_pad(c1, settings.out_dir, q2_hist + "_bkg_subtracted.pdf");
    //bg->Draw("SAME");
    TH1* data_subtracted = subtract_bkg(h_data, bg);

    TF1* sgnl = new TF1("sgnl", normed_gauss_binwidth, signal_min, signal_max, 4);
    sgnl->SetParameters(10, 960, 5, 10);
    sgnl->SetParNames("Number events", "mean", "sigma", "bin width");
    sgnl->SetParLimits(0, 0, 1000);
    sgnl->SetParLimits(1, 958, 962);
    sgnl->SetParLimits(2, 1, 15);
    //sgnl->FixParameter(1, 960);  // obtained by summing up all IMee bins and fitting the resulting signal peak
    sgnl->FixParameter(3, h_data->GetBinWidth(h_data->FindFixBin(958)));
    data_subtracted->Fit(sgnl, Form("%sLB", default_fit_options));  // B bounds (use given parameter limits)

    if (!settings.out_dir.empty()) {
        TF1* tot = new TF1("tot", total, fit_min, fit_max, 8);
        tot->SetParameters(bg->GetParameter(0), bg->GetParameter(1), bg->GetParameter(2), bg->GetParameter(3),
                           sgnl->GetParameter(0), sgnl->GetParameter(1), sgnl->GetParameter(2), sgnl->GetParameter(3));
        tot->SetNpx(1000);
        tot->SetLineColor(kRed+1);
        c1->cd();
        h_data->Draw();  // draw this one new for an empty histogram
        h_data->GetXaxis()->SetRangeUser(combined_plot_xRange.Start(), combined_plot_xRange.Stop());
        tot->Draw("SAME");
        bg->SetLineColor(kGreen+2);
        bg->Draw("SAME");
        c1->Modified();
        c1->Print((settings.out_dir / fs::path(q2_hist + "_combined.pdf")).c_str());
    }

    c2->cd();
    data_subtracted->Draw();
    sgnl->SetLineColor(kRed+1);
    sgnl->Draw("SAME");
    LOG(INFO) << "Fitted number of signal events for bin " << q2_hist << ": " << sgnl->GetParameter(0) << " +/- " << sgnl->GetParError(0);
    imee_fits.emplace_back(sgnl->GetParameter(0), sgnl->GetParError(0));
    // only show the fit information on the second canvas
    auto ps = dynamic_cast<TPaveStats*>(c2->GetPrimitive("stats"));
    ps->SetOptFit(111);
    c2->Update();
    save_pad(c2, settings.out_dir, q2_hist + "_signal_fit.pdf");

    }  // end loop imee_bins



    // default numbers, TODO: obtain automatically from provided MC file
    vector<int> N_true = {1796465, 648784, 443837, 413906, 321715, 266582, 268359, 269363, 273838, 367069, 336559, 301793, 203500, 37855};
    vector<int> N_rec = {44647, 82368, 82391, 76605, 60388, 49237, 51209, 46409, 42275, 53206, 48704, 43092, 27068, 6236};

    for (auto i : imee_bins) {
        double eff = N_rec.at(i);
        eff /= N_true.at(i);
        double eff_err = eff * (1 - eff) / N_true.at(i);

        double corr = imee_fits.at(i).first / eff;
        double corr_err = 1./eff * sqrt(std_ext::sqr(imee_fits.at(i).second) + std_ext::sqr(imee_fits.at(i).first) * std_ext::sqr(eff_err));
        signal_fits.emplace_back(corr, corr_err);

        LOG(INFO) << "Efficiency corrected signal events bin " << hist_names.at(i) << ": " << corr << " +/- " << corr_err;
    }
    if (el::Loggers::verboseLevel()) {
        auto sum_pair = [] (double sum, const pair<double, double> n) {
            return sum + n.first;
        };
        LOG(INFO) << "Total number of fitted signal events: " << accumulate(signal_fits.begin(), signal_fits.end(), 0., sum_pair);
        double N_rec_total = accumulate(N_rec.begin(), N_rec.end(), 0.), N_true_total = accumulate(N_true.begin(), N_true.end(), 0.);
        LOG(INFO) << "Overall analysis efficiency: " << N_rec_total / N_true_total << "%   (" << N_rec_total << "/" << N_true_total << ")";
    }
}


template <typename T>
struct TCLAPInterval : interval<T> {
    using interval<T>::interval;
    using ValueCategory = TCLAP::ValueLike;
};

int main(int argc, char** argv) {
    SetupLogger();

    // parse command line
    TCLAP::CmdLine cmd("EtapDalitz_fit", ' ', "0.1");
    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"int");
    auto cmd_batchmode = cmd.add<TCLAP::MultiSwitchArg>("b","batch","Run in batch mode (no ROOT shell afterwards)",false);
    auto cmd_debug = cmd.add<TCLAP::MultiSwitchArg>("","debug","Enable debug mode",false);

    auto cmd_ref = cmd.add<TCLAP::MultiSwitchArg>("r","reference","Run Reference Channel Analysis", false);
    auto cmd_ref_only = cmd.add<TCLAP::MultiSwitchArg>("","ref-only","Only Reference Channel Analysis", false);

    auto cmd_input = cmd.add<TCLAP::ValueArg<string>>("i","input","ROOT input file",true,"","rootfile");
    auto cmd_mcinput = cmd.add<TCLAP::ValueArg<string>>("m","mcinput","Input for MC histograms",false,"","rootfile");
    auto cmd_output = cmd.add<TCLAP::ValueArg<string>>("o","output","Output file",false,"","filename");
    auto cmd_out_dir = cmd.add<TCLAP::ValueArg<string>>("d","directory","Output directory for images/PDFs",false,"","directory");
    TCLAP::ValuesConstraintExtra<decltype(ExpConfig::Setup::GetNames())> allowedsetupnames(ExpConfig::Setup::GetNames());
    auto cmd_setup  = cmd.add<TCLAP::ValueArg<string>>("s","setup","Choose setup by name",false,"", &allowedsetupnames);

    auto cmd_imee_bins = cmd.add<TCLAP::ValueArg<string>>("","bins","Comma-separated ranges of q2 bins, no spaces, e.g. 2-6,9",false,
                                                          string("0-")+to_string(q2_params_t::bin_widths.size()-1),"bins");
    auto cmd_EPTrange = cmd.add<TCLAP::ValueArg<string>>("c","EPTrange","EPT channel range for reference fits, e.g. 0-40 or 0-10,35-40",
                                                         false,"0-40","channels");
    auto cmd_rebin = cmd.add<TCLAP::ValueArg<int>>("","rebin","Number of bins to rebin for some to-be-fitted histograms",false,0,"int");

    cmd.parse(argc, argv);

    constexpr int etap_threshold_eptCh = 40;

    const bool ref = cmd_ref->isSet();
    const bool ref_only = cmd_ref_only->isSet();

    // retrieve settings_t singleton and set values (defaults defined via TCLAP)
    settings_t& settings = settings_t::get();
    settings.rebin = cmd_rebin->getValue();
    if (settings.rebin < 0) {
        LOG(WARNING) << "Provided rebin value is negative! rebin will be ignored.";
        settings.rebin = 0;
    }
    if (settings.rebin)
        LOG(INFO) << "Some of the to-be-fitted histograms will be rebinned combining " << settings.rebin << " bins";

    settings.out_dir = cmd_out_dir->getValue();
    if (!settings.out_dir.empty() && !fs::exists(settings.out_dir)) {
        LOG(WARNING) << "The specified output directory doesn't exist! It will be created.";
        if (!fs::create_directories(settings.out_dir))
            LOG(FATAL) << "Failed to create output directory.";
    } else if (fs::exists(settings.out_dir) && !fs::is_directory(settings.out_dir))
        LOG(FATAL) << "The specified output directory is not a directory!";
    if (!settings.out_dir.empty())
        LOG(INFO) << "Images/PDFs will be saved in the following output directory: " << settings.out_dir.string();

    const bool debug = settings.debug = cmd_debug->isSet();
    // verbosity management
    if (debug)
        el::Loggers::setVerboseLevel(1);  // if debug is chosen, show at least some debug output from the logger if no other verbosity level is provided
    if (cmd_verbose->isSet()) {
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());
    }
    // silence RooFit's messenger
    if (cmd_verbose->getValue() < 3)
        RooMsgService::instance().setGlobalKillBelow(debug ? WARNING : ERROR);
    if (!debug)
        RooMsgService::instance().setSilentMode(true);
    // silence unwanted ROOT info messages like Canvas::Print file has been created etc.
    if (!debug)
        gErrorIgnoreLevel = kWarning;

    // do some tests in the beginning to make sure all functions work as expected
    if (debug) {
        test_path_building();
        cout << "\nCall cut extraction method\n" << endl;
        print_extracted_cuts(cmd_input->getValue());
    }

    if (cmd_setup->isSet())
        ExpConfig::Setup::SetByName(cmd_setup->getValue());
    else {
        LOG(WARNING) << "No Setup specified, use \"Setup_2014_07_EPT_Prod\" as default fallback";
        ExpConfig::Setup::SetByName("Setup_2014_07_EPT_Prod");
    }

    WrapTFileInput input(cmd_input->getValue());
    WrapTFileInput mcinput;
    if (cmd_mcinput->isSet())
        mcinput.OpenFile(cmd_mcinput->getValue());

    auto taggChRange = convert_piecewise_interval_type<int>(
                progs::tools::parse_cmdline_ranges(std_ext::tokenize_string(cmd_EPTrange->getValue(), ",")), true);
    if (cmd_EPTrange->isSet())
        LOG(WARNING) << "Using non-default Tagger channel range, may not yield correct results (debugging purposes)";
    // tagger channel range of interest: 0 - 40 (where 40 contains the eta' threshold)
    // make sure the constructed channel list matches these conditions
    if (taggChRange.back() > etap_threshold_eptCh) {
        LOG(WARNING) << "Highest EPT channel " << taggChRange.back() << " provided is below the eta' threshold! "
                     << "All channels above " << etap_threshold_eptCh << " will be skipped";
        auto pos = find(taggChRange.begin(), taggChRange.end(), etap_threshold_eptCh);
        auto elements = distance(pos, taggChRange.end());
        taggChRange.erase(pos, taggChRange.end());
        VLOG(1) << "Removed " << elements << " elements from EPT channel range";
        VLOG(3) << "New range is now: " << taggChRange;
    }
    // reverse the Tagger channel range to loop over them later in increasing beam energy steps
    reverse(taggChRange.begin(), taggChRange.end());

    const auto q2_bins = convert_piecewise_interval(progs::tools::parse_cmdline_ranges(std_ext::tokenize_string(cmd_imee_bins->getValue(), ",")));
    if (q2_bins.back() >= q2_params_t::bin_widths.size())
        LOG(FATAL) << "Provided range for q2 bins goes up to " << q2_bins.back() << ", but largest bin index is " << q2_params_t::bin_widths.size()-1;
    if (debug)
        cout << "parsed the following bins: " << q2_bins << endl;

    // create TRint as RooFit internally creates functions/histograms,
    // prevents this stupid gStyle=0 related error, sigh...
    argc = 0;  // prevent TRint to parse any cmdline
    // IMPORTANT! Create TRint on the heap to prevent ROOT from segfaulting when closing the ROOT shell
    auto app = new TRint("EtapDalitz_fit", &argc, argv, nullptr, 0, true);
    if (cmd_batchmode->isSet())
        gROOT->SetBatch(true);

    // set signal handler after starting TRint, otherwise it will be overwritten with ROOT handlers
    signal(SIGINT, [] (int) {
        LOG(WARNING) << "Processing interrupted";
        interrupt = true;
    });

    unique_ptr<WrapTFileOutput> masterFile;
    if (cmd_output->isSet()) {
        // cd into masterFile upon creation
        masterFile = std_ext::make_unique<WrapTFileOutput>(cmd_output->getValue(), true);
    }

    // globally applied style settings
    // change line scaling for PDF output, default is three (which results in ugly thick lines)
    gStyle->SetLineScalePS(2);
    // change title size, especially important for RooFit since I couldn't figure out how else to change it
    gStyle->SetTitleFontSize(.07f);


    pair<double, double> ref_fit;
    vector<pair<double, double>> signal_fits = {};

    if (ref || ref_only) {
        TH1* ref;
        input.GetObject("EtapDalitz_plot_Ref/KinFitProb > 0.01/PID E cut < 0.3 MeV/h/Data/etapIM_kinfitted", ref);

        /* code using double gaussian
        // fit ARGUS model with double gauss for reference
        // --- Observable ---
        RooRealVar mes("IM","IM_{#gamma#gamma}", 840, 1020, "MeV");

        // --- Parameters ---
        RooRealVar sigmean1("sigmean1","#eta' mass", 958., 950., 965.);
        RooRealVar sigwidth1("sigwidth1","#eta' width", 2., .1, 30.);
        RooRealVar sigmean2("sigmean2","#eta' mass", 959., 950., 965.);
        RooRealVar sigwidth2("sigwidth2","#eta' width", 3., .1, 30.);

        // --- Build Double Gaussian PDF ---
        RooGaussian gauss1("gauss1", "gauss1 PDF", mes, sigmean1, sigwidth1);
        RooGaussian gauss2("gauss2", "gauss2 PDF", mes, sigmean2, sigwidth2);
        RooRealVar g1frac("g1frac","fraction of gauss1",.4,0.,1.);
        RooAddPdf signalModel("doubleGauss","g1+g2", RooArgList(gauss1,gauss2), g1frac);

        // --- Build Argus background PDF ---
        RooRealVar argpar("argpar","argus shape parameter",-5.,-25.,5.);
        RooArgusBG background("background","Argus PDF",mes,RooConst(1000),argpar);

        // --- Construct signal+background PDF ---
        RooRealVar nsig("nsig","#signal events",1000,0.,1e5);
        RooRealVar nbkg("nbkg","#background events",1000,0.,1e5);
        RooAddPdf model("model","g+a",RooArgList(signalModel,background),RooArgList(nsig,nbkg));

        RooDataHist h_roo_data("h_roo_data","dataset",mes,ref);
        */
        // fit ARGUS model with CB function for reference
        // --- Observable ---
        RooRealVar mes("IM","IM_{#gamma#gamma}", 840, 1020, "MeV");

        // --- Build CB Function PDF ---
        RooRealVar cb_x0("cb_x0", "expectation value", 958, 950, 970);
        RooRealVar cb_sigma("cb_sigma", "standard deviation", 2, .01, 20);
        RooRealVar cb_alpha("cb_alpha", "transition gauss to power function", 1, .1, 20);
        RooRealVar cb_n("cb_n", "parameter power function", 1, .1, 5);
        RooCBShape signal("signal", "CB Function", mes, cb_x0, cb_sigma, cb_alpha, cb_n);

        // --- Build Argus background PDF ---
        RooRealVar argpar("argpar","argus shape parameter",-5.,-25.,5.);
        RooArgusBG background("background","Argus PDF",mes,RooConst(1000),argpar);

        // --- Construct signal+background PDF ---
        RooRealVar nsig("nsig","#signal events",1000,0.,1e5);
        RooRealVar nbkg("nbkg","#background events",1000,0.,1e5);
        RooAddPdf model("model","g+a",RooArgList(signal,background),RooArgList(nsig,nbkg));

        RooDataHist h_roo_data("h_roo_data","dataset",mes,ref);

        // --- Perform extended ML fit of composite PDF to data ---
        model.fitTo(h_roo_data, Extended(), SumW2Error(kTRUE));

        // --- Plot toy data and composite PDF overlaid ---
        RooPlot* mesframe = mes.frame();
        h_roo_data.plotOn(mesframe);
        model.plotOn(mesframe);
        model.plotOn(mesframe, Components(background), LineStyle(ELineStyle::kDashed));

        mesframe->Draw();
        gPad->Modified();
        gPad->Update();

        reference_fit(input, "KinFitProb > 0.01/PID E cut < 0.3 MeV", taggChRange, mcinput, ref_fit);
    } else
        ref_fit = {6.00768e+06, 59855.2};  // use default values for the fit result if reference_fit wasn't used

    if (!ref_only && !interrupt)
        signal_fit(input, {}, q2_bins, mcinput, signal_fits);


    // run TRint
    if (!cmd_batchmode->isSet()) {
        if (!std_ext::system::isInteractive()) {
            LOG(INFO) << "No TTY attached. Not starting ROOT shell.";
        }
        else {
            if (masterFile)
                LOG(INFO) << "Close ROOT properly to write data to disk.";

            app->Run(kTRUE);  // really important to return...
            if (masterFile)
                LOG(INFO) << "Writing output file...";
            masterFile = nullptr;  // and to destroy the master WrapTFile before TRint is destroyed
            // call this before application tear down
            gROOT->EndOfProcessCleanups();
            // do not delete app, otherwise ROOT might segfault
        }
    }

    return EXIT_SUCCESS;
}
