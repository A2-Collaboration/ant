/**
  * @file EtapDalitz_plot.cc
  * @brief Plotting tool for TTrees from the ant::analysis::physics::etaprime::etaprime_dalitz physics class
  */

#include "base/Logger.h"

#include "analysis/plot/CutTree.h"
#include "analysis/utils/particle_tools.h"
#include "analysis/physics/etaprime/etaprime_dalitz.h"

#include "base/CmdLine.h"
#include "base/interval.h"
#include "base/printable.h"
#include "base/WrapTFile.h"
#include "base/std_ext/string.h"
#include "base/std_ext/system.h"
#include "base/std_ext/math.h"
#include "base/vec/vec2.h"
#include "base/ProgressCounter.h"

#include "expconfig/ExpConfig.h"

#include "tree/TAntHeader.h"

#include <list>
#include <vector>
#include <algorithm>

#include "TSystem.h"
#include "TRint.h"
#include "TStyle.h"

using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::plot;
using namespace std;

volatile sig_atomic_t interrupt = false;
static double binScale = 1.;


template<typename Hist_t>
struct MCTrue_Splitter : cuttree::StackedHists_t<Hist_t> {

    // Hist_t should have that type defined
    using Fill_t = typename Hist_t::Fill_t;

    const decltype(physics::EtapDalitz::makeChannels()) channels;

    constexpr static Color_t bkg_color = kGray+1;

    MCTrue_Splitter(const HistogramFactory& histFac,
                    const cuttree::TreeInfo_t& treeInfo) :
        cuttree::StackedHists_t<Hist_t>(histFac, treeInfo),
        channels(physics::EtapDalitz::makeChannels())
    {
        using histstyle::Mod_t;

        this->GetHist(0, "Data",   Mod_t::MakeDataPoints(kBlack));
        this->GetHist(1, "Signal", Mod_t::MakeLine(kRed, 2));
        // mctrue is never >= 3 (and < 9) in tree, use this to sum up all MC and all bkg MC
        // see also Fill()
        this->GetHist(2, "Sum_MC", Mod_t::MakeLine(kBlack, 1));
        this->GetHist(3, "Bkg_MC", Mod_t::MakeFill(bkg_color, -1));
    }

    void Fill(const Fill_t& f) {

        const unsigned mctrue = unsigned(f.Tree.channel);

        //if (mctrue >= 10 && !(mctrue == 10 || mctrue == 12))
        //    return;

        using histstyle::Mod_t;

        auto get_bkg_name = [] (const unsigned mctrue) {
            const auto entry = physics::EtapDalitz::reaction_channels.channels.find(int(mctrue));

            if (entry != physics::EtapDalitz::reaction_channels.channels.end())
                return entry->second.name;

            return string("Unknown Decay");
        };

        auto get_color = [] (const unsigned mctrue) -> short {
            const auto entry = physics::EtapDalitz::reaction_channels.channels.find(int(mctrue));

            if (entry != physics::EtapDalitz::reaction_channels.channels.end())
                return entry->second.color;

            return histstyle::color_t::Get(mctrue-10);
        };

        const Hist_t& hist = mctrue < 10 ? this->GetHist(mctrue)
                                         : this->GetHist(mctrue,
                                                         get_bkg_name(mctrue),
                                                         Mod_t::MakeLine(get_color(mctrue), 1, bkg_color)
                                                         );

        hist.Fill(f);

        // handle MC_all and MC_bkg
        if (mctrue > 0) {
            this->GetHist(2).Fill(f);
            if (mctrue >= 10)
                this->GetHist(3).Fill(f);
        }
    }
};

bool Contains(const interval<double>& i, const std::vector<double>& d)
{
    for (const auto& v : d)
        if (i.Contains(v))
            return true;

    return false;
}

double max(const std::vector<double>& data)
{
    return *max_element(data.cbegin(), data.cend());
}

double im_ee(vector<double> vetoE, vector<TLorentzVector> photons)
{
    auto find_leptons = [] (vector<double> vec) {
        vector<size_t> p(vec.size());
        std::iota(p.begin(), p.end(), 0);
        std::sort(p.begin(), p.end(),
                  [&] (size_t i, size_t j) {
            return vec[i] > vec[j];
        });
        return p;
    };

    auto leptons = find_leptons(vetoE);

    return (photons.at(leptons[0]) + photons.at(leptons[1])).M();
}


// define the structs containing the histograms and the cuts
struct Hist_t {

    using Tree_t = physics::EtapDalitz::Tree_t;

    struct Fill_t {
        const Tree_t& Tree;

        Fill_t(const Tree_t& t) : Tree(t) {}

        double TaggW() const {
            return Tree.TaggW;
        }
    };

    template <typename Hist>
    using fillfunc_t = std::function<void(Hist*, const Fill_t&)>;

    template <typename Hist>
    struct HistFiller_t {
        fillfunc_t<Hist> func;
        Hist* h;
        HistFiller_t(Hist* hist, fillfunc_t<Hist> f): func(f), h(hist) {}
        void Fill(const Fill_t& data) const {
            func(h, data);
        }
    };

    template <typename Hist>
    struct HistMgr : std::list<HistFiller_t<Hist>> {

        using list<HistFiller_t<Hist>>::list;

        void Fill(const Fill_t& data) const {
            for (auto& h : *this)
                h.Fill(data);
        }
    };

    static BinSettings Bins(const unsigned bins, const double min, const double max) {
        return BinSettings(unsigned(bins*binScale), min, max);
    }

    HistMgr<TH1D> h1;
    HistMgr<TH2D> h2;

    const BinSettings Ebins    = Bins(1200, 0, 1200);

    const BinSettings Chi2Bins = BinSettings(250, 0, 25);
    const BinSettings probbins = BinSettings(250, 0,  1);

    const BinSettings IMbins   = Bins(1200,   0, 1200);
    const BinSettings MMbins   = Bins(1200, 400, 1600);

    const BinSettings TaggChBins = BinSettings(47);

    const BinSettings TaggTime   = BinSettings(240, -30, 30);
    const BinSettings CoplBins   = Bins(300, 0, 30);

    const BinSettings zVertex    = Bins(100, -15, 15);
    const BinSettings free_vz    = Bins(400, -40, 40);

    HistogramFactory HistFac;

    void AddTH1(const string &title, const string &xlabel, const string &ylabel,
                const BinSettings &bins, const string &name, fillfunc_t<TH1D> f)
    {
        h1.emplace_back(HistFiller_t<TH1D>(
                            HistFac.makeTH1D(title, xlabel, ylabel, bins, name), f));
    }

    void AddTH2(const string &title, const string &xlabel, const string &ylabel,
                const BinSettings &xbins, const BinSettings& ybins,
                const string &name, fillfunc_t<TH2D> f)
    {
        h2.emplace_back(HistFiller_t<TH2D>(
                            HistFac.makeTH2D(title, xlabel, ylabel, xbins, ybins, name), f));
    }

    Hist_t(const HistogramFactory& hf, cuttree::TreeInfo_t): HistFac(hf) {

        AddTH1("KinFitChi2", "#chi^{2}", "#", Chi2Bins, "KinFitChi2",
               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.kinfit_chi2, f.TaggW());
        });

        AddTH1("TreeFitChi2", "#chi^{2}", "#", Chi2Bins, "TreeFitChi2",
               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.treefit_chi2, f.TaggW());
        });

        AddTH1("KinFitProb", "probability", "#", probbins, "KinFitProb",
               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.kinfit_probability, f.TaggW());
        });

        AddTH1("TreeFitProb", "probability", "#", probbins, "TreeFitProb",
               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.treefit_probability, f.TaggW());
        });

        AddTH1("3 photon IM", "3#gamma IM [MeV]", "#", IMbins, "etapIM",
               [] (TH1D* h, const Fill_t& f) {
            h->Fill(f.Tree.etap().M(), f.TaggW());
        });

        AddTH1("3 photon IM kinfitted",  "3#gamma IM fit [MeV]", "#", IMbins, "etapIM_kinfitted",
               [] (TH1D* h, const Fill_t& f) {
            h->Fill(f.Tree.etap_kinfit().M(), f.TaggW());
        });

        AddTH1("Missing Mass", "MM [MeV]", "", MMbins, "mm",
               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.mm().M(), f.TaggW());
        });

        AddTH1("Z Vertex Kinfit", "z [cm]", "#", zVertex, "v_z_kinfit",
               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.kinfit_ZVertex, f.TaggW());
        });

        AddTH1("Z Vertex free Z Kinfit", "z [cm]", "#", free_vz, "v_z_kinfit_freeZ",
               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.kinfit_freeZ_ZVertex, f.TaggW());
        });

        AddTH1("Z Vertex Treefit", "z [cm]", "#", zVertex, "v_z_treefit",
               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.treefit_ZVertex, f.TaggW());
        });

//        AddTH2("IM(e+e-) vs. IM(e+e-g)", "IM(e+e-g) [MeV]", "IM(e+e-) [MeV]", BinSettings(600, 0, 1200), BinSettings(500, 0, 1000), "IM2d",
//               [] (TH2D* h, const Fill_t& f) {
//            h->Fill(f.Tree.etap().M(), im_ee(f.Tree.photons_vetoE(), f.Tree.photons()), f.TaggW());
//        });

//        AddTH2("IM(e+e-) vs. IM(e+e-g) prompt", "IM(e+e-g) [MeV]", "IM(e+e-) [MeV]", BinSettings(600, 0, 1200), BinSettings(500, 0, 1000), "IM2d_prompt",
//               [] (TH2D* h, const Fill_t& f) {
//            if (f.TaggW() > 0)
//                h->Fill(f.Tree.etap().M(), im_ee(f.Tree.photons_vetoE(), f.Tree.photons()));
//        });

//        AddTH2("IM(e+e-) vs. IM(e+e-g) random", "IM(e+e-g) [MeV]", "IM(e+e-) [MeV]", BinSettings(600, 0, 1200), BinSettings(500, 0, 1000), "IM2d_random",
//               [] (TH2D* h, const Fill_t& f) {
//            if (f.TaggW() < 0)
//                h->Fill(f.Tree.etap().M(), im_ee(f.Tree.photons_vetoE(), f.Tree.photons()));
//        });

//        AddTH2("IM(e+e-) vs. IM(e+e-g) fit", "IM(e+e-g) [MeV]", "IM(e+e-) [MeV]", BinSettings(600, 0, 1200), BinSettings(500, 0, 1000), "IM2d_fit",
//               [] (TH2D* h, const Fill_t& f) {
//            h->Fill(f.Tree.etap_kinfit().M(), im_ee(f.Tree.photons_vetoE(), f.Tree.photons_kinfitted()), f.TaggW());
//        });
        });

//        AddTH2("IM(e+e-) vs. IM(e+e-g) fit prompt", "IM(e+e-g) [MeV]", "IM(e+e-) [MeV]", BinSettings(600, 0, 1200), BinSettings(500, 0, 1000), "IM2d_fit_prompt",
//               [] (TH2D* h, const Fill_t& f) {
//            if (f.TaggW() > 0)
//                h->Fill(f.Tree.etap_kinfit().M(), im_ee(f.Tree.photons_vetoE(), f.Tree.photons_kinfitted()));
//        });

//        AddTH2("IM(e+e-) vs. IM(e+e-g) fit random", "IM(e+e-g) [MeV]", "IM(e+e-) [MeV]", BinSettings(600, 0, 1200), BinSettings(500, 0, 1000), "IM2d_fit_random",
//               [] (TH2D* h, const Fill_t& f) {
//            if (f.TaggW() < 0)
//                h->Fill(f.Tree.etap_kinfit().M(), im_ee(f.Tree.photons_vetoE(), f.Tree.photons_kinfitted()));
//        });

//        AddTH2("Cluster Size vs. Energy", "Energy [MeV]", "Cluster Size", Ebins, BinSettings(50), "clusterSize_E",
//               [] (TH2D* h, const Fill_t& f) {
//            for (unsigned i = 0; i < f.Tree.photons().size(); i++)
//                h->Fill(f.Tree.photons().at(i).Energy(), f.Tree.photons_clusterSize().at(i), f.TaggW());
//        });

        AddTH1("Effective Cluster Radius", "R", "#", BinSettings(500, 0, 50), "clusterRadius",
               [] (TH1D* h, const Fill_t& f) {
            for (unsigned i = 0; i < f.Tree.photons().size(); i++)
                h->Fill(f.Tree.photons_effect_radius().at(i), f.TaggW());
        });

        AddTH1("Lateral Moment", "L", "#", BinSettings(200, 0, 1), "lateralMoment",
               [] (TH1D* h, const Fill_t& f) {
            for (unsigned i = 0; i < f.Tree.photons().size(); i++)
                h->Fill(f.Tree.photons_lat_moment().at(i), f.TaggW());
        });

        AddTH2("Effective Cluster Radius vs. Energy", "Energy [MeV]", "R", Bins(300, 0, 1200), BinSettings(200, 0, 50), "clusterRadius_E",
               [] (TH2D* h, const Fill_t& f) {
            for (unsigned i = 0; i < f.Tree.photons().size(); i++)
                h->Fill(f.Tree.photons().at(i).Energy(), f.Tree.photons_effect_radius().at(i), f.TaggW());
        });

        AddTH2("Lateral Moment vs. Energy", "Energy [MeV]", "L", Bins(300, 0, 1200), BinSettings(100, 0, 1), "lateralMoment_E",
               [] (TH2D* h, const Fill_t& f) {
            for (unsigned i = 0; i < f.Tree.photons().size(); i++)
                h->Fill(f.Tree.photons().at(i).Energy(), f.Tree.photons_lat_moment().at(i), f.TaggW());
        });

//        AddTH2("Lateral Moment vs. Effective Cluster Radius", "R", "L", BinSettings(500, 0, 50), BinSettings(200, 0, 1), "lateralMoment_clusterRadius",
//               [] (TH2D* h, const Fill_t& f) {
//            for (unsigned i = 0; i < f.Tree.photons().size(); i++)
//                h->Fill(f.Tree.photons_effect_radius().at(i), f.Tree.photons_lat_moment().at(i), f.TaggW());
//        });

    }

    void Fill(const Fill_t& f) const
    {
        h1.Fill(f);
        h2.Fill(f);
    }

    std::vector<TH1*> GetHists() const
    {
        vector<TH1*> v;
        v.reserve(h1.size()+h2.size());
        for (auto& e : h1)
            v.emplace_back(e.h);
        for (auto& e: h2)
            v.emplace_back(e.h);
        return v;
    }

    // Sig and Ref channel share some cuts...
    static cuttree::Cuts_t<Fill_t> GetCuts()
    {

        using cuttree::MultiCut_t;

        cuttree::Cuts_t<Fill_t> cuts;

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"KinFitProb>0.001", [] (const Fill_t& f) { return f.Tree.kinfit_probability > .001; }},
                              {"KinFitProb>0.02", [] (const Fill_t& f) { return f.Tree.kinfit_probability > .02; }},
                              {"KinFitProb>0.05", [] (const Fill_t& f) { return f.Tree.kinfit_probability > .05; }}
                             });

        auto antiPi0Cut = [] (const Fill_t& f, const double low = 102., const double high = 170.) {
            const interval<double> pion_cut(low, high);
            TLorentzVector pi0;
            const std::vector<std::array<size_t, 2>> pi0_combs = {{0, 2}, {1, 2}};

            auto sort_energies = [] (vector<double> vec) {
                vector<size_t> p(vec.size());
                std::iota(p.begin(), p.end(), 0);
                std::sort(p.begin(), p.end(),
                          [&] (size_t i, size_t j) {
                    return vec[i] > vec[j];
                });
                return p;
            };

            const auto photons = f.Tree.photons();
            const auto sorted = sort_energies(f.Tree.photons_vetoE());

            for (const auto pi0_comb : pi0_combs) {
                pi0 = TLorentzVector(0., 0., 0., 0.);

                for (const auto idx : pi0_comb)
                    pi0 += TParticle(ParticleTypeDatabase::Photon, photons.at(sorted.at(idx)));

                // check anti pi^0 cut
                if (pion_cut.Contains(pi0.M()))
                    return false;
            }

            return true;
        };

        auto distinctPIDCut = [] (const Fill_t& f) {
            const auto photons = f.Tree.photons();
            const auto vetos = f.Tree.photons_vetoE();
            const auto channels = f.Tree.photons_vetoChannel();

            const unsigned size = photons.size();
            unsigned idx[size];
            for (unsigned i = 0; i < size; i++)
                idx[i] = i;

            // sort idx according to the photons' veto energies
            sort(idx, idx+size, [vetos] (unsigned i, unsigned j) { return vetos[i] > vetos[j]; });

            return channels.at(idx[0]) != channels.at(idx[1]);
        };

        auto freeZ_vertexCut = [] (const Fill_t& f, const double high = 6.) {
            return f.Tree.kinfit_freeZ_ZVertex < high;
        };

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"anti pi0", antiPi0Cut}
                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"distinct PID elements", distinctPIDCut}
                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"free vz cut", freeZ_vertexCut}
                          });

        return cuts;
    }

};


int main(int argc, char** argv)
{
    SetupLogger();

    TCLAP::CmdLine cmd("plot", ' ', "0.1");
    auto cmd_input = cmd.add<TCLAP::ValueArg<string>>("i", "input", "Input file", true, "", "input");
    auto cmd_batchmode = cmd.add<TCLAP::MultiSwitchArg>("b", "batch", "Run in batch mode (no ROOT shell afterwards)", false);
    auto cmd_maxevents = cmd.add<TCLAP::MultiArg<int>>("m", "maxevents", "Process only max events", false, "maxevents");
    auto cmd_output = cmd.add<TCLAP::ValueArg<string>>("o", "output", "Output file", false, "", "filename");

    auto cmd_tree = cmd.add<TCLAP::ValueArg<string>>("", "tree", "Tree name", false, "T1", "treename");
    auto cmd_pres = cmd.add<TCLAP::SwitchArg>("p", "", "Presentation Mode", false);
    auto cmd_binscale = cmd.add<TCLAP::ValueArg<double>>("B", "bin-scale", "Bin Scaling", false, 1.0, "bins");

    cmd.parse(argc, argv);


    // open the TRint app as early as possible to prevent ROOT to create a new one automatically
    // which will cause problems because of bad ROOT internal data / pointer handling, might cause segfaults
    argc = 0;  // prevent TRint to parse any cmdline
    TRint app("EtapDalitz_plot", &argc, argv, nullptr, 0, true);


    // set signal handler after starting TRint, otherwise it will be overwritten with ROOT handlers
    signal(SIGINT, [] (int) {
        LOG(WARNING) << "Processing interrupted";
        interrupt = true;
    });


    // general styling settings for transparence
//    gStyle->SetFillColor(0);
//    gStyle->SetFillStyle(0);
//    gStyle->SetFrameFillColor(0);
//    gStyle->SetFrameFillStyle(0);

    if (cmd_pres->isSet()) {
        gStyle->SetLabelSize(.05f, "XYZ");
        gStyle->SetTextSize(.05f);
        gStyle->SetCanvasBorderSize(0);
    }

    if (cmd_binscale->isSet()) {
        binScale = cmd_binscale->getValue();
    }


    WrapTFileInput input;
    string setup_name;
    try {

        input.OpenFile(cmd_input->getValue());

        auto header = input.GetSharedClone<TAntHeader>("AntHeader");

        if (!header) {
            LOG(WARNING) << "No TAntHeader found in " << cmd_input->getValue();
            return 1;
        }

        setup_name = header->SetupName;

    } catch (const std::runtime_error& e) {
        LOG(ERROR) << "Can't open " << cmd_input->getValue() << " " << e.what();
    }

    if (setup_name.empty())
        setup_name = "Setup_2014_07_EPT_Prod";

    ExpConfig::Setup::Get(setup_name);


    auto link_branches = [&input] (const string treename, WrapTTree* wraptree, long long expected_entries) {
        TTree* t;
        if (!input.GetObject(treename,t))
            throw runtime_error("Cannot find tree " + treename + " in input file");
        if (expected_entries >= 0 && t->GetEntries() != expected_entries)
            throw runtime_error("Tree " + treename + " does not have entries == " + to_string(expected_entries));
        if (wraptree->Matches(t, false)) {
            wraptree->LinkBranches(t);
            return true;
        }
        return false;
    };


    Hist_t::Tree_t tree;

    if (!link_branches("EtapDalitz/tree", addressof(tree), -1)) {
        LOG(ERROR) << "Cannot link branches of tree";
        return 1;
    }

    const auto entries = tree.Tree->GetEntries();

    unique_ptr<WrapTFileOutput> masterFile;
    if (cmd_output->isSet()) {
        masterFile = std_ext::make_unique<WrapTFileOutput>(cmd_output->getValue(),
                                                           WrapTFileOutput::mode_t::recreate,
                                                           true); // cd into masterFile upon creation
    }


    HistogramFactory HistFac("EtapDalitz");

    const auto& sanitized_treename = std_ext::replace_str(cmd_tree->getValue(),"/","_");

    auto signal_hists = cuttree::Make<MCTrue_Splitter<Hist_t>>(HistFac,
                                                               sanitized_treename,
                                                               Hist_t::GetCuts()
                                                               );

    LOG(INFO) << "Tree entries = " << entries;

    auto max_entries = entries;
    if (cmd_maxevents->isSet() && cmd_maxevents->getValue().back() < entries) {
        max_entries = cmd_maxevents->getValue().back();
        LOG(INFO) << "Running until " << max_entries;
    }

    long long entry = 0;
    double last_percent = 0;
    ProgressCounter::Interval = 3;
    ProgressCounter progress(
                [&entry, entries, &last_percent] (std::chrono::duration<double> elapsed) {
        const double percent = 100.*entry/entries;
        const double speed = (percent - last_percent)/elapsed.count();
        LOG(INFO) << setw(2) << setprecision(4) << "Processed " << percent << " %, ETA: " << ProgressCounter::TimeToStr((100-percent)/speed);
        last_percent = percent;
    });

    while (entry < max_entries) {
        if (interrupt)
            break;

        tree.Tree->GetEntry(entry++);
        /*if (tree.MCtrue && tree.CBSumE < 550.)
            continue;
        const double taggTime = tree.TaggT - tree.CBAvgTime;
        //std::cout << "tagg time: " << taggTime << std::endl;
        if (tree.TaggW > 0 && (taggTime > 1.5 || taggTime < -2.))  // tighten prompt peak -> TaggW not right anymore!!!
            continue;*/
        cuttree::Fill<MCTrue_Splitter<Hist_t>>(signal_hists, {tree});

        //entry++;

        ProgressCounter::Tick();
    }

    LOG(INFO) << "Analyzed " << entry << " events, speed "
              << entry/progress.GetTotalSecs() << " event/s";

    if (!cmd_batchmode->isSet()) {
        if (!std_ext::system::isInteractive())
            LOG(INFO) << "No TTY attached. Not starting ROOT shell.";
        else {
            if (masterFile)
                LOG(INFO) << "Stopped running, but close ROOT properly to write data to disk.";

            app.Run(kTRUE); // really important to return...
            if (masterFile)
                LOG(INFO) << "Writing output file...";
            masterFile = nullptr;   // and to destroy the master WrapTFile before TRint is destroyed
        }
    }

    return EXIT_SUCCESS;
}
