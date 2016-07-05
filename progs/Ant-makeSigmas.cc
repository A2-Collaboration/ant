

#include "analysis/utils/PullsWriter.h"
#include "analysis/plot/HistogramFactories.h"

#include "base/std_ext/system.h"
#include "base/Logger.h"
#include "base/CmdLine.h"
#include "base/WrapTFile.h"
#include "base/ProgressCounter.h"
#include "base/std_ext/string.h"

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

    for(int z = 0; z < bins.Bins(); ++z) {
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

void ClearHistogram(TH2D* hist) {
    for(int x=0; x<=hist->GetNbinsX(); ++x) {
        for(int y=0; y<=hist->GetNbinsY(); ++y) {
            hist->SetBinContent(x,y, std_ext::NaN);
        }
    }
}

void FitSlicesZ(const TH3D* hist, HistogramFactory& hf_, const bool do_fit=false) {

    HistogramFactory hf(formatter() << hist->GetName() << "_FitZ", hf_);

    const auto xbins = getBins(hist->GetXaxis());
    const auto ybins = getBins(hist->GetYaxis());
    const auto zbins = getBins(hist->GetZaxis());

    vector<TH2D*> parmhists(3, nullptr);

    for(size_t p=0; p<parmhists.size(); ++p) {
        parmhists.at(p) = hf.makeTH2D(
                              formatter() << hist->GetTitle() << ": Parameter " << p,
                              hist->GetXaxis()->GetTitle(),
                              hist->GetYaxis()->GetTitle(),
                              xbins,
                              ybins,
                              formatter() << hist->GetName() << "_z_Parm" << p );
    }

    TH2D* h_chi2 = hf.makeTH2D(
                       formatter() << hist->GetTitle() << ":Chi2",
                       hist->GetXaxis()->GetTitle(),
                       hist->GetYaxis()->GetTitle(),
                       xbins,
                       ybins,
                       formatter() << hist->GetName() << "_z_Chi2" );

    TH2D* h_RMS  = hf.makeTH2D(
                       formatter() << hist->GetTitle() << ": RMS",
                       hist->GetXaxis()->GetTitle(),
                       hist->GetYaxis()->GetTitle(),
                       xbins,
                       ybins,
                       formatter() << hist->GetName() << "_z_RMS" );
    ClearHistogram(h_RMS);

    TH2D* h_Mean  = hf.makeTH2D(
                       formatter() << hist->GetTitle() << ": Mean",
                       hist->GetXaxis()->GetTitle(),
                       hist->GetYaxis()->GetTitle(),
                       xbins,
                       ybins,
                       formatter() << hist->GetName() << "_z_Mean" );
    ClearHistogram(h_Mean);

    TCanvas* c = new TCanvas();
    const string c_title = formatter() << hist->GetTitle() << " Fits";
    c->SetTitle(c_title.c_str());
    c->Divide(xbins.Bins(), ybins.Bins());


    for(int x=0; x<xbins.Bins(); ++x) {
        for(int y=0; y<ybins.Bins(); ++y) {

            auto pad = c->cd(1+x+(ybins.Bins()-y-1)*xbins.Bins());

            //pad->SetMargin(0.01,0.01,0.01,0.01);

            TH1D* slice = projectZ(hist, x+1, y+1, hf);

            slice->Draw();

            if(slice->Integral() > 500.0) {
                const auto rms = slice->GetRMS();
                const auto mean = slice->GetMean();

                h_RMS->SetBinContent(x+1,y+1, rms);
                h_Mean->SetBinContent(x+1,y+1, mean);

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

                    h_chi2->SetBinContent(x+1,y+1, chi2);

                    if(chi2 > .0) {

                        if(tf1_gaus->GetNpar() != parmhists.size())
                            throw std::runtime_error("Wrong number pf parameters");

                        for(size_t p=0; p<parmhists.size(); ++p) {
                            auto h =  parmhists.at(p);
                            h->SetBinContent(x+1, y+1, tf1_gaus->GetParameter(p));
                            h->SetBinError(x+1,y+1,tf1_gaus->GetParError(p));
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
}



int main( int argc, char** argv )
{
    SetupLogger();

    signal(SIGINT, [] (int) {
        LOG(INFO) << ">>> Interrupted";
        interrupt = true;
    });

    TCLAP::CmdLine cmd("Ant-makeSigmas", ' ', "0.1");

    auto cmd_input       = cmd.add<TCLAP::ValueArg<string>>("i","input","pull trees",true,"","rootfile");
    auto cmd_output      = cmd.add<TCLAP::ValueArg<string>>("o","","sigma hists",false,"","rootfile");
    auto cmd_batchmode   = cmd.add<TCLAP::MultiSwitchArg>("b","batch","Run in batch mode (no ROOT shell afterwards)",false);
    auto cmd_maxevents   = cmd.add<TCLAP::MultiArg<int>>("m","maxevents","Process only max events",false,"maxevents");
    auto cmd_tree        = cmd.add<TCLAP::ValueArg<string>>("t","tree","Treename",false,"JustPi0/m2Pi0/pulls_photon_cb","treename");
    auto cmd_fitprob_cut = cmd.add<TCLAP::ValueArg<double>>("","fitprob_cut","Min. required Fit Probability",false,0.01,"probability");

    cmd.parse(argc, argv);

    const auto fitprob_cut = cmd_fitprob_cut->getValue();

    WrapTFileInput input(cmd_input->getValue());

    TTree* tree;
    if(!input.GetObject(cmd_tree->getValue(), tree)) {
        LOG(ERROR) << "Cannot find tree in " << cmd_input->getValue();
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
            s.bins_E = {10, 0, 800};
        }
        if(std_ext::string_ends_with(treename, "pulls_photon_taps")) {
            s.name = "sigma_photon_taps";
            s.bins_cosTheta = bins_cosTheta_taps;
            s.bins_E = {10, 0, 800};
        }
        if(std_ext::string_ends_with(treename, "pulls_proton_cb")) {
            s.name = "sigma_proton_cb";
            s.bins_cosTheta = bins_cosTheta_cb;
            s.bins_E = {5, 0, 400};
        }
        if(std_ext::string_ends_with(treename, "pulls_proton_taps")) {
            s.name = "sigma_proton_taps";
            s.bins_cosTheta = bins_cosTheta_taps;
            s.bins_E = {5, 0, 400};
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

            FitSlicesZ(h_pullsE, HistFac);
            FitSlicesZ(h_pullsTheta, HistFac);
            FitSlicesZ(h_pullsPhi, HistFac);

            FitSlicesZ(h_sigmasE, HistFac);
            FitSlicesZ(h_sigmasTheta, HistFac);
            FitSlicesZ(h_sigmasPhi, HistFac);

            if(masterFile)
                LOG(INFO) << "Stopped running, but close ROOT properly to write data to disk.";

            app.Run(kTRUE); // really important to return...
            if(masterFile)
                LOG(INFO) << "Writing output file...";
            masterFile = nullptr;   // and to destroy the master WrapTFile before TRint is destroyed
        }
    }

}
