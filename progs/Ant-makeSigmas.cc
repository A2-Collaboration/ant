

#include "analysis/utils/PullsWriter.h"
#include "analysis/plot/HistogramFactories.h"

#include "base/std_ext/system.h"
#include "base/Logger.h"
#include "base/CmdLine.h"
#include "base/WrapTFile.h"
#include "base/ProgressCounter.h"
#include "base/std_ext/string.h"
#include "base/Array2D.h"

#include "analysis/plot/root_draw.h"

#include "TTree.h"
#include "TRint.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TCanvas.h"

using namespace ant;
using namespace std;
using namespace ant::analysis;
using namespace ant::std_ext;
static volatile bool interrupt = false;


bool IsBinValid(const TH2D* hist, const int x, const int y) {
    return (x>0 && x <=hist->GetNbinsX()) && (y>0 && y <=hist->GetNbinsY());
}

const std::vector<std::pair<int,int>> neighbors = {{+1,0},{-1,0},{0,+1},{0,-1}};

double getNeighborAverage(const TH2D* hist, const int x, const int y) {

    double sum = {};
    unsigned n = 0;

    for(const auto& d : neighbors) {
        const auto bx = x + d.first;
        const auto by = y + d.second;

        if(IsBinValid(hist,bx,by)) {
            const auto v = hist->GetBinContent(bx,by);
            if(std::isfinite(v)) {
                sum += v;
                ++n;
            }
        }
    }

    return n>0 ? sum / n : 0.0;
}

unsigned getNeighborCount(const TH2D* hist, const int x, const int y) {

    unsigned n = 0;

    for(const auto& d : neighbors) {
        const auto bx = x + d.first;
        const auto by = y + d.second;

        const auto valid = IsBinValid(hist,bx,by);
        if( valid && isfinite(hist->GetBinContent(bx,by))) {
            ++n;
        }
    }

    return n;
}

void fillNeighborAverages(TH2D* hist) {


    unsigned neighbors=0;

    do {
        neighbors=0;
        int p_x =0;
        int p_y =0;

        for(int x=1; x<=hist->GetNbinsX(); ++x) {
            for(int y=1; y<=hist->GetNbinsY(); ++y) {

                if(std::isnan(hist->GetBinContent(x,y))) {
                    const auto n = getNeighborCount(hist, x,y);
                    if(n>neighbors) {
                        neighbors = n;
                        p_x = x;
                        p_y = y;
                    }
                }
            }
        }

        // if updatable bin found
        if(neighbors > 0) {
            const auto a = getNeighborAverage(hist,p_x,p_y);
            hist->SetBinContent(p_x,p_y, a);
        }
    } while(neighbors != 0);
}



BinSettings getBins(const TAxis* axis) {
    return BinSettings(axis->GetNbins(), axis->GetXmin(), axis->GetXmax());
}

/**
 * @brief projectZ
 *        code after TH3::FitSlicesZ()
 * @param hist
 * @param x
 * @param y
 * @param hf
 * @return
 *
 */

TH1D* projectZ(const TH3D* hist, const int x, const int y, HistogramFactory& hf) {

    const string name = formatter() << hist->GetName() << "_z_" << x << "_" << y;

    const auto axis = hist->GetZaxis();
    const auto bins = getBins(axis);

    auto h = hf.makeTH1D(name.c_str(), hist->GetZaxis()->GetTitle(), "", bins, name.c_str());
    h->Reset();

    for(int z = 0; z < int(bins.Bins()); ++z) {
        const auto v = hist->GetBinContent(x, y, z);
        if(v != .0) {
            h->Fill(axis->GetBinCenter(z), v);
            h->SetBinError(z, hist->GetBinError(x,y,z));
        }
    }

    return h;
}

vec2 maximum(const TH1D* hist) {
    return vec2( hist->GetBinCenter(hist->GetMaximumBin()), hist->GetMaximum() );
}

void ClearHistogram(TH2D* hist, const double v=0.0) {
    for(int x=1; x<=hist->GetNbinsX(); ++x) {
        for(int y=1; y<=hist->GetNbinsY(); ++y) {
            hist->SetBinContent(x,y, v);
        }
    }
}

interval<double> getZRange(const TH2& hist) {
    return {hist.GetMinimum(), hist.GetMaximum()};
}

void MakeSameZRange(std::vector<TH2*> hists) {

    if(hists.empty())
        return;


    auto b1 = getZRange(*hists.front());

    for(auto it=next(hists.begin()); it!=hists.end(); ++it) {
        b1.Extend(getZRange(**it));
    }

    for(auto h : hists) {
        h->SetMinimum(b1.Start());
        h->SetMaximum(b1.Stop());
    }

}

struct FitSlices1DHists {
    TH2D* Mean    = nullptr;
    TH2D* RMS     = nullptr;
    TH2D* Entries = nullptr;

    TH2D* chi2dof = nullptr;


    std::vector<TH2D*> parameter = {};
};

FitSlices1DHists FitSlicesZ(const TH3D* hist, HistogramFactory& hf_, const bool do_fit=false, const string& title="", const double min_intragral=1000.0) {

    HistogramFactory hf(formatter() << hist->GetName() << "_FitZ", hf_);

    const auto xbins = getBins(hist->GetXaxis());
    const auto ybins = getBins(hist->GetYaxis());
    const auto zbins = getBins(hist->GetZaxis());

    FitSlices1DHists result;

    result.parameter.reserve(3);

    for(size_t p=0; p<result.parameter.size(); ++p) {
        result.parameter.at(p) = hf.makeTH2D(
                              formatter() << hist->GetTitle() << ": Parameter " << p,
                              hist->GetXaxis()->GetTitle(),
                              hist->GetYaxis()->GetTitle(),
                              xbins,
                              ybins,
                              formatter() << hist->GetName() << "_z_Parm" << p );
    }

    result.chi2dof  = hf.makeTH2D(
                       formatter() << hist->GetTitle() << ":Chi2",
                       hist->GetXaxis()->GetTitle(),
                       hist->GetYaxis()->GetTitle(),
                       xbins,
                       ybins,
                       formatter() << hist->GetName() << "_z_Chi2" );
    result.chi2dof->SetStats(false);

    result.RMS  = hf.makeTH2D(
                       formatter() << hist->GetTitle() << ": RMS",
                       hist->GetXaxis()->GetTitle(),
                       hist->GetYaxis()->GetTitle(),
                       xbins,
                       ybins,
                       formatter() << hist->GetName() << "_z_RMS" );
    ClearHistogram(result.RMS, std_ext::NaN);
    result.RMS->SetStats(false);

    result.Mean  = hf.makeTH2D(
                       formatter() << hist->GetTitle() << ": Mean",
                       hist->GetXaxis()->GetTitle(),
                       hist->GetYaxis()->GetTitle(),
                       xbins,
                       ybins,
                       formatter() << hist->GetName() << "_z_Mean" );
    ClearHistogram(result.Mean, std_ext::NaN);
    result.Mean->SetStats(false);

    result.Entries  = hf.makeTH2D(
                       formatter() << hist->GetTitle() << ": Entries",
                       hist->GetXaxis()->GetTitle(),
                       hist->GetYaxis()->GetTitle(),
                       xbins,
                       ybins,
                       formatter() << hist->GetName() << "_z_Entries" );
    result.Entries->SetStats(false);

    TCanvas* c = new TCanvas();
    const string c_title = formatter() << title << ": " << hist->GetTitle() << " Fits";
    c->SetTitle(c_title.c_str());
    c->Divide(int(xbins.Bins()), int(ybins.Bins()));


    for(int x=0; x < int(xbins.Bins()); ++x) {
        for(int y=0; y < int(ybins.Bins()); ++y) {

            const int padno = int(1+x+(ybins.Bins()-y-1)*xbins.Bins());
            auto pad = c->cd(padno);

            //pad->SetMargin(0.01,0.01,0.01,0.01);

            TH1D* slice = projectZ(hist, x+1, y+1, hf);

            slice->Draw();

            const auto integral = slice->Integral();
            result.Entries->SetBinContent(x+1,y+1, integral);

            if(integral > min_intragral) {
                const auto rms = slice->GetRMS();
                const auto mean = slice->GetMean();

                result.RMS->SetBinContent(x+1,y+1, rms);
                result.Mean->SetBinContent(x+1,y+1, mean);

                if(do_fit) {

                    auto tf1_gaus = new TF1("gaus","gaus");
                    tf1_gaus->SetNpx(500);

                    const auto max = maximum(slice);

                    tf1_gaus->SetParameter(0, max.y);
                    tf1_gaus->SetParLimits(0, .1, 2.*max.y);

                    tf1_gaus->SetParameter(1, max.x);

                    tf1_gaus->SetParameter(2, rms);
                    tf1_gaus->SetParLimits(2, 0.0, zbins.Length());

                    slice->Fit(tf1_gaus, "BQ");

                    const auto chi2 = tf1_gaus->GetChisquare() / tf1_gaus->GetNDF();

                    result.chi2dof->SetBinContent(x+1,y+1, chi2);

                    if(chi2 > .0) {

                        if(size_t(tf1_gaus->GetNpar()) != result.parameter.size())
                            throw std::runtime_error("Wrong number pf parameters");

                        for(size_t p=0; p<result.parameter.size(); ++p) {
                            auto h =  result.parameter.at(p);
                            h->SetBinContent(x+1, y+1, tf1_gaus->GetParameter(int(p)));
                            h->SetBinError  (x+1, y+1, tf1_gaus->GetParError (int(p)));
                        }
                    } else {
                        pad->SetFillColor(kYellow);
                    }

                }
            } else {
                pad->SetFillColor(kGray);
            }

        }
    }

    return result;
}

struct NewSigamsResult_t {
    TH2D* newSigmas = nullptr;
    TH2D* oldSigmas = nullptr;
    TH2D* pulls     = nullptr;
};

NewSigamsResult_t makeNewSigmas(const TH3D* pulls, const TH3D* sigmas, HistogramFactory& hf, const string& output_name, const string& treename, const double integral_cut) {
    const string newTitle = formatter() << "new " << sigmas->GetTitle();

    auto pull_values  = FitSlicesZ(pulls,  hf, false, treename, integral_cut);
    auto sigma_values = FitSlicesZ(sigmas, hf, false, treename, integral_cut);

    NewSigamsResult_t result;

    result.oldSigmas = sigma_values.Mean;
    result.pulls     = pull_values.RMS;

    result.newSigmas = hf.clone(result.oldSigmas, output_name);
    assert(result.newSigmas);

    result.newSigmas->Multiply(result.pulls);

    //fillNeighborAverages(result.newSigmas);
    {
        auto wrapper = Array2D_TH2D(result.newSigmas);
        FloodFillAverages::fillNeighborAverages(wrapper);
    }

    result.newSigmas->SetTitle(newTitle.c_str());

    MakeSameZRange( {result.newSigmas,result.oldSigmas} );

    return result;
}



int main( int argc, char** argv )
{
    SetupLogger();

    signal(SIGINT, [] (int) {
        LOG(INFO) << ">>> Interrupted";
        interrupt = true;
    });

    TCLAP::CmdLine cmd("Ant-makeSigmas", ' ', "0.1");

    auto cmd_input        = cmd.add<TCLAP::ValueArg<string>>("i","input",      "pull trees",                                     true,  "", "rootfile");
    auto cmd_output       = cmd.add<TCLAP::ValueArg<string>>("o","",           "sigma hists",                                    false, "", "rootfile");
    auto cmd_batchmode    = cmd.add<TCLAP::MultiSwitchArg>  ("b","batch",      "Run in batch mode (no ROOT shell afterwards)",   false);
    auto cmd_maxevents    = cmd.add<TCLAP::MultiArg<int>>   ("m","maxevents",  "Process only max events",                        false, "maxevents");
    auto cmd_tree         = cmd.add<TCLAP::ValueArg<string>>("t","tree",       "Treename",false,"JustPi0/m2Pi0/pulls_photon_cb","treename");
    auto cmd_fitprob_cut  = cmd.add<TCLAP::ValueArg<double>>("", "fitprob_cut","Min. required Fit Probability",                  false, 0.01,"probability");
    auto cmd_integral_cut = cmd.add<TCLAP::ValueArg<double>>("", "integral_cut","Min. required integral in Bins",                false, 100.0,"integral");

    cmd.parse(argc, argv);

    const auto fitprob_cut  = cmd_fitprob_cut->getValue();
    const auto integral_cut = cmd_integral_cut->getValue();

    WrapTFileInput input(cmd_input->getValue());

    TTree* tree;
    if(!input.GetObject(cmd_tree->getValue(), tree)) {
        LOG(ERROR) << "Cannot find tree " << cmd_tree->getValue() << " in " << cmd_input->getValue();
        exit(EXIT_FAILURE);
    }

    utils::PullsWriter::PullTree_t pulltree;
    if(!pulltree.Matches(tree)) {
        LOG(ERROR) << "Given tree is not a PullTree_t";
        exit(EXIT_FAILURE);
    }
    pulltree.LinkBranches(tree);
    auto entries = pulltree.Tree->GetEntries();

    unique_ptr<WrapTFileOutput> masterFile;
    if(cmd_output->isSet()) {
        masterFile = std_ext::make_unique<WrapTFileOutput>(cmd_output->getValue(),
                                                    WrapTFileOutput::mode_t::recreate,
                                                     true); // cd into masterFile upon creation
    }

    struct hist_settings_t {
        string name;
        BinSettings bins_cosTheta{0};
        BinSettings bins_E{0};
    };

    auto get_hist_settings = [] (const std::string& treename) {
        hist_settings_t s;
        BinSettings bins_cosTheta_cb{15, std::cos(std_ext::degree_to_radian(160.0)), std::cos(std_ext::degree_to_radian(20.0))};
        BinSettings bins_cosTheta_taps{10, std::cos(std_ext::degree_to_radian(25.0)), std::cos(std_ext::degree_to_radian(0.0))};
        if(std_ext::string_ends_with(treename, "pulls_photon_cb")) {
            s.name = "sigma_photon_cb";
            s.bins_cosTheta = bins_cosTheta_cb;
            s.bins_E = {10, 0, 1000};
        }
        if(std_ext::string_ends_with(treename, "pulls_photon_taps")) {
            s.name = "sigma_photon_taps";
            s.bins_cosTheta = bins_cosTheta_taps;
            s.bins_E = {10, 0, 1000};
        }
        if(std_ext::string_ends_with(treename, "pulls_proton_cb")) {
            s.name = "sigma_proton_cb";
            s.bins_cosTheta = bins_cosTheta_cb;
            s.bins_E = {1, 0, 1000};
        }
        if(std_ext::string_ends_with(treename, "pulls_proton_taps")) {
            s.name = "sigma_proton_taps";
            s.bins_cosTheta = bins_cosTheta_taps;
            s.bins_E = {1, 0, 1000};
        }
        return s;
    };

    const hist_settings_t& hist_settings = get_hist_settings(cmd_tree->getValue());
    if(hist_settings.name.empty()) {
        LOG(ERROR) << "Could not identify pulltree name";
        exit(EXIT_FAILURE);
    }
    HistogramFactory HistFac(hist_settings.name);

    BinSettings bins_pulls(50,-5,5);

    auto h_pullsE = HistFac.makeTH3D("Pulls Ek","cos #theta","Ek","Pulls Ek",
                                     hist_settings.bins_cosTheta,
                                     hist_settings.bins_E,
                                     bins_pulls,
                                     "h_pullsE"
                                     );

    auto h_pullsTheta = HistFac.makeTH3D("Pulls Theta","cos #theta","Ek","Pulls #theta",
                                     hist_settings.bins_cosTheta,
                                     hist_settings.bins_E,
                                     bins_pulls,
                                     "h_pullsTheta"
                                     );

    auto h_pullsPhi = HistFac.makeTH3D("Pulls Phi","cos #theta","Ek","Pulls #phi",
                                     hist_settings.bins_cosTheta,
                                     hist_settings.bins_E,
                                     bins_pulls,
                                     "h_pullsPhi"
                                     );

    auto h_sigmasE = HistFac.makeTH3D("Sigmas Ek","cos #theta","Ek","Sigmas Ek",
                                     hist_settings.bins_cosTheta,
                                     hist_settings.bins_E,
                                     BinSettings(50,0.0,50.0),
                                     "h_sigmasE"
                                     );

    auto h_sigmasTheta = HistFac.makeTH3D("Sigmas Theta","cos #theta","Ek","Sigmas #theta",
                                     hist_settings.bins_cosTheta,
                                     hist_settings.bins_E,
                                     BinSettings(50,0.0,degree_to_radian(10.0)),
                                     "h_sigmasTheta"
                                     );

    auto h_sigmasPhi = HistFac.makeTH3D("Sigmas Phi","cos #theta","Ek","Sigmas #phi",
                                     hist_settings.bins_cosTheta,
                                     hist_settings.bins_E,
                                     BinSettings(50,0.0,degree_to_radian(10.0)),
                                     "h_sigmasPhi"
                                     );


    LOG(INFO) << "Tree entries=" << entries;
    auto max_entries = entries;
    if(cmd_maxevents->isSet() && cmd_maxevents->getValue().back()<entries) {
        max_entries = cmd_maxevents->getValue().back();
        LOG(INFO) << "Running until " << max_entries;
    }

    long long entry = 0;
    ProgressCounter::Interval = 3;
    ProgressCounter progress(
                [&entry, entries] (std::chrono::duration<double>) {
        LOG(INFO) << "Processed " << 100.0*entry/entries << " %";
    });

    LOG(INFO) << " Min. required probability = " << fitprob_cut;

    for(entry=0;entry<max_entries;entry++) {
        if(interrupt)
            break;

        progress.Tick();
        pulltree.Tree->GetEntry(entry);

        if(pulltree.FitProb > fitprob_cut ) {

            h_pullsE->Fill(     cos(pulltree.Theta), pulltree.E, pulltree.PullE,      pulltree.TaggW);
            h_pullsTheta->Fill( cos(pulltree.Theta), pulltree.E, pulltree.PullTheta,  pulltree.TaggW);
            h_pullsPhi->Fill(   cos(pulltree.Theta), pulltree.E, pulltree.PullPhi,    pulltree.TaggW);

            h_sigmasE->Fill(    cos(pulltree.Theta), pulltree.E, pulltree.SigmaE,     pulltree.TaggW);
            h_sigmasTheta->Fill(cos(pulltree.Theta), pulltree.E, pulltree.SigmaTheta, pulltree.TaggW);
            h_sigmasPhi->Fill(  cos(pulltree.Theta), pulltree.E, pulltree.SigmaPhi,   pulltree.TaggW);

        }

    }

    if(!cmd_batchmode->isSet()) {
        if(!std_ext::system::isInteractive()) {
            LOG(INFO) << "No TTY attached. Not starting ROOT shell.";
        }
        else {
            argc=1; // prevent TRint to parse any cmdline except prog name
            TRint app("Ant-makeSigmas",&argc,argv,nullptr,0,true);

            auto new_E     = makeNewSigmas(h_pullsE,     h_sigmasE,     HistFac, "sigma_E",    cmd_tree->getValue(), integral_cut);
            auto new_Theta = makeNewSigmas(h_pullsTheta, h_sigmasTheta, HistFac, "sigma_Theta",cmd_tree->getValue(), integral_cut);
            auto new_Phi   = makeNewSigmas(h_pullsPhi,   h_sigmasPhi,   HistFac, "sigma_Phi",  cmd_tree->getValue(), integral_cut);

            canvas summary(cmd_tree->getValue());
            summary << drawoption("colz");

            for( auto r : {new_E, new_Theta, new_Phi}) {
                summary << r.newSigmas << r.oldSigmas << r.pulls;
            }
            summary << endc;

            if(masterFile)
                LOG(INFO) << "Stopped running, but close ROOT properly to write data to disk.";

            app.Run(kTRUE); // really important to return...
            if(masterFile)
                LOG(INFO) << "Writing output file...";
            masterFile = nullptr;   // and to destroy the master WrapTFile before TRint is destroyed
        }
    }

}
