

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
#include "TF1.h"
#include "TCanvas.h"

using namespace ant;
using namespace std;
using namespace ant::analysis;
using namespace ant::std_ext;
static volatile bool interrupt = false;

BinSettings getBins(const TAxis* axis) {
    return BinSettings(axis->GetNbins(), axis->GetXmin(), axis->GetXmax());
}

TH1D* projectZ(const TH3D* hist, const int x, const int y, HistogramFactory& hf) {

    const string name = formatter() << hist->GetName() << "_z_" << x << "_" << y;

    const auto bins = getBins(hist->GetZaxis());

    auto h = hf.makeTH1D(name.c_str(), hist->GetZaxis()->GetTitle(), "", bins, name.c_str());

    for(int z = 0; z < bins.Bins(); ++z) {
        h->SetBinContent(z, hist->GetBinContent(x,y, z));
        h->SetBinError(z, hist->GetBinError(x,y,z));
    }

    return h;
}

void FitSlicesZ(const TH3D* hist, HistogramFactory& hf_) {

    HistogramFactory hf(formatter() << hist->GetName() << "_FitZ", hf_);

    const int nx = hist->GetNbinsX();
    const int ny = hist->GetNbinsY();

    TCanvas* c = new TCanvas();
    c->Divide(nx,ny);

    int pad_number = 1;
    for(int x=0; x<nx; ++x) {
        for(int y=0; y<ny; ++y) {

            c->cd(pad_number++);

            TH1D* slice = projectZ(hist, x, y, hf);

            //slice->Draw();

            auto tf1_gaus = new TF1("gaus","gaus");
            tf1_gaus->SetParameter(0, slice->GetMaximum());
            tf1_gaus->SetParameter(1, slice->GetBinCenter(slice->GetMaximumBin()));
            tf1_gaus->SetParameter(2, slice->GetRMS());

            slice->Fit(tf1_gaus);

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

    auto cmd_input = cmd.add<TCLAP::ValueArg<string>>("i","input","pull trees",true,"","rootfile");
    auto cmd_output = cmd.add<TCLAP::ValueArg<string>>("o","","sigma hists",false,"","rootfile");
    auto cmd_batchmode = cmd.add<TCLAP::MultiSwitchArg>("b","batch","Run in batch mode (no ROOT shell afterwards)",false);
    auto cmd_maxevents = cmd.add<TCLAP::MultiArg<int>>("m","maxevents","Process only max events",false,"maxevents");
    auto cmd_tree = cmd.add<TCLAP::ValueArg<string>>("t","tree","Treename",false,"JustPi0/m2Pi0/pulls_photon_cb","treename");


    cmd.parse(argc, argv);

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

    for(entry=0;entry<max_entries;entry++) {
        if(interrupt)
            break;

        pulltree.Tree->GetEntry(entry);

        if(pulltree.FitProb > .1 ) {

            h_pullsE->Fill(cos(pulltree.Theta), pulltree.E, pulltree.PullE, pulltree.TaggW);
            h_pullsTheta->Fill(cos(pulltree.Theta), pulltree.E, pulltree.PullTheta, pulltree.TaggW);
            h_pullsPhi->Fill(cos(pulltree.Theta), pulltree.E, pulltree.PullPhi, pulltree.TaggW);

        }

    }

    // important to call FitSlicesZ with already created function,
    // otherwise gROOT global is used which forces the creation of default TApplication
    // leading to weird side effects with other globals such as gStyle when instantiating TRint below
   // auto tf1_gaus = new TF1("gaus","gaus");

    //h_pullsE->FitSlicesZ(tf1_gaus);
//    h_pullsTheta->FitSlicesZ(tf1_gaus);
//    h_pullsPhi->FitSlicesZ(tf1_gaus);



    if(!cmd_batchmode->isSet()) {
        if(!std_ext::system::isInteractive()) {
            LOG(INFO) << "No TTY attached. Not starting ROOT shell.";
        }
        else {
            argc=1; // prevent TRint to parse any cmdline except prog name
            TRint app("Ant-makeSigmas",&argc,argv,nullptr,0,true);

            FitSlicesZ(h_pullsTheta, HistFac);

            if(masterFile)
                LOG(INFO) << "Stopped running, but close ROOT properly to write data to disk.";

            app.Run(kTRUE); // really important to return...
            if(masterFile)
                LOG(INFO) << "Writing output file...";
            masterFile = nullptr;   // and to destroy the master WrapTFile before TRint is destroyed
        }
    }

}
