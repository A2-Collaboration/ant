#include "base/Logger.h"

#include "analysis/plot/CutTree.h"
#include "analysis/physics/etaprime/etaprime_omega_gamma.h"

#include "base/CmdLine.h"
#include "base/interval.h"
#include "base/printable.h"
#include "base/WrapTFile.h"
#include "base/std_ext/string.h"
#include "base/std_ext/system.h"

#include "TSystem.h"
#include "TRint.h"

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::plot;
using namespace std;

volatile bool interrupt = false;

template<typename Hist_t>
struct MCTrue_Splitter : cuttree::StackedHists_t<Hist_t> {

    // Hist_t should have that type defined
    using Fill_t = typename Hist_t::Fill_t;

    MCTrue_Splitter(const HistogramFactory& histFac) : cuttree::StackedHists_t<Hist_t>(histFac) {
        using histstyle::Mod_t;
        this->GetHist(0, "Data", Mod_t::MakeDataPoints(kBlack));
        this->GetHist(1, "Sig",  Mod_t::MakeLine(kRed, 2.0));
        this->GetHist(2, "Ref",  Mod_t::MakeLine(kRed, 2.0));
        // mctrue is never >=3 (and <9) in tree, use this to sum up all MC and all bkg MC
        // see also Fill()
        this->GetHist(3, "Sum_MC", Mod_t::MakeLine(kBlack, 2.0));
        this->GetHist(4, "Bkg_MC", Mod_t::MakeFill(kGray+2));
    }

    void Fill(const Fill_t& f) {

        const unsigned mctrue = f.Common.MCTrue;

        auto get_bkg_name = [] (unsigned mctrue) {
            const string& name = mctrue>=10 ?
                                     physics::EtapOmegaG::ptreeBackgrounds[mctrue-10].Name
                                 : "Other";
            return "Bkg_"+name;
        };

        using histstyle::Mod_t;
        const Hist_t& hist = mctrue<9 ? this->GetHist(mctrue) :
                                        this->GetHist(mctrue,
                                                      get_bkg_name(mctrue),
                                                      Mod_t::MakeLine(histstyle::color_t::Get(mctrue-9), 1, kGray+2)
                                                      );

        hist.Fill(f);

        // handle MC_all and MC_bkg
        if(mctrue>0) {
            this->GetHist(3).Fill(f);
            if(mctrue >= 9)
                this->GetHist(4).Fill(f);
        }
    }
};

// define the structs containing the histograms
// and the cuts. for simple branch variables, that could
// be combined...

struct CommonHist_t {
    using Tree_t = physics::EtapOmegaG::TreeCommon;
    struct Fill_t {
        const Tree_t& Common;
        Fill_t(const Tree_t& common) : Common(common) {}
        double TaggW() const {
            return Common.TaggW;
        }
    };



    const BinSettings bins_FitProb{200, 0, 0.3};
    TH1D* h_KinFitProb;
    TH1D* h_CBSumE;
    TH1D* h_CBSumVetoE;
    TH1D* h_PIDSumE;


    CommonHist_t(HistogramFactory HistFac) {
        h_KinFitProb = HistFac.makeTH1D("KinFitProb","p","",bins_FitProb,"h_KinFitProb");
        h_CBSumE = HistFac.makeTH1D("CB Sum E","E / MeV","",BinSettings(100,500,1600),"h_CBSumE");
        h_CBSumVetoE = HistFac.makeTH1D("CB Veto Sum E","E / MeV","",BinSettings(50,0,10),"h_CBSumVetoE");
        h_PIDSumE = HistFac.makeTH1D("PID Sum E","E / MeV","",BinSettings(50,0,10),"h_PIDSumE");
    }
    void Fill(const Fill_t& f) const {
        h_CBSumE->Fill(f.Common.CBSumE, f.TaggW());
        h_CBSumVetoE->Fill(f.Common.CBSumVetoE, f.TaggW());
        h_PIDSumE->Fill(f.Common.PIDSumE, f.TaggW());
    }
    std::vector<TH1*> GetHists() const {
        return {h_KinFitProb, h_CBSumE, h_CBSumVetoE, h_PIDSumE};
    }

    // Sig and Ref channel share some cuts...
    static cuttree::Cuts_t<Fill_t> GetCuts() {
        using cuttree::MultiCut_t;
        cuttree::Cuts_t<Fill_t> cuts;
        cuts.emplace_back(MultiCut_t<Fill_t>{
                              // Use non-null PID cuts only when PID calibrated...
                              //{"CBSumVeto=0", [] (const Fill_t& f) { return f.Common.CBSumVetoE==0; } },
                              {"CBSumVeto<0.25", [] (const Fill_t& f) { return f.Common.CBSumVetoE<0.25; } },
                              //{"PIDSumE=0", [] (const Fill_t& f) { return f.Common.PIDSumE==0; } },
                              {"PIDSumE<1", [] (const Fill_t& f) { return f.Common.PIDSumE<1; } },
                          });
//        cuts.emplace_back(MultiCut_t<Fill_t>{
//                                 {"KinFitProb>0", [] (const Fill_t& f) { return f.Common.KinFitProb>0; } },
//                                 {"KinFitProb>0.01", [] (const Fill_t& f) { return f.Common.KinFitProb>0.01; } },
//                             });
        return cuts;
    }


};

struct SigHist_t : CommonHist_t {
    using SharedTree_t = physics::EtapOmegaG::Sig_t::SharedTree_t;
    using Tree_t = physics::EtapOmegaG::Sig_t::Fit_t::BaseTree_t;

    struct Fill_t : CommonHist_t::Fill_t {
        const SharedTree_t& Shared;
        const Tree_t& Tree;
        Fill_t(const CommonHist_t::Tree_t& common,
               const SharedTree_t& shared,
               const Tree_t& tree) :
            CommonHist_t::Fill_t(common),
            Shared(shared),
            Tree(tree)
        {}
    };

    TH2D* h_IM_gg_gg;     // Goldhaber plot
    TH1D* h_IM_4g;        // EtaPrime IM
    TH1D* h_IM_gg;        // EtaPrime IM

    TH1D* h_AntiPi0FitProb;
    TH1D* h_AntiEtaFitProb;
    TH1D* h_TreeFitProb;

    const BinSettings bins_IM_Etap {100, 800,1050};
    const BinSettings bins_IM_Omega{100, 550, 950};
    const BinSettings bins_ClusterShape{30,0,1};

    SigHist_t(HistogramFactory HistFac) : CommonHist_t(HistFac) {
        BinSettings bins_goldhaber(200, 0, 900);
        const string axislabel_goldhaber("2#gamma IM / MeV");

        h_IM_gg_gg = HistFac.makeTH2D("IM 2#gamma-2#gamma",
                                    axislabel_goldhaber, axislabel_goldhaber,
                                    bins_goldhaber, bins_goldhaber,
                                    "h_IM_gg_gg"
                                    );
        h_IM_4g = HistFac.makeTH1D("#eta' IM", "IM(#pi^{0}#gamma#gamma) / MeV","",bins_IM_Etap,"h_IM_4g");
        h_IM_gg = HistFac.makeTH1D("#eta' IM", "IM(#gamma#gamma) / MeV","",BinSettings(200,80,700),"h_IM_gg");
        h_AntiPi0FitProb = HistFac.makeTH1D("AntiPi0FitProb", "p","",bins_FitProb,"h_AntiPi0FitProb");
        h_AntiEtaFitProb = HistFac.makeTH1D("AntiEtaFitProb", "p","",bins_FitProb,"h_AntiEtaFitProb");
        h_TreeFitProb = HistFac.makeTH1D("TreeFitProb", "p","",bins_FitProb,"h_TreeFitProb");

    }

    void Fill(const Fill_t& f) const {
        CommonHist_t::Fill(f);
        const SharedTree_t& s = f.Shared;
        const Tree_t& tree = f.Tree;

        for(unsigned i=0;i<s.gg_gg1().size();i++) {
            h_IM_gg_gg->Fill(s.gg_gg1()[i], s.gg_gg2()[i], f.TaggW());
            h_IM_gg_gg->Fill(s.gg_gg2()[i], s.gg_gg1()[i], f.TaggW());
        }
        h_AntiPi0FitProb->Fill(s.AntiPi0FitProb, f.TaggW());
        h_AntiEtaFitProb->Fill(s.AntiEtaFitProb, f.TaggW());

        h_IM_4g->Fill(tree.IM_Pi0gg, f.TaggW());
        h_IM_gg->Fill(tree.IM_gg, f.TaggW());
        h_TreeFitProb->Fill(tree.TreeFitProb, f.TaggW());

    }

    std::vector<TH1*> GetHists() const {
        auto hists = CommonHist_t::GetHists();
        hists.insert(hists.end(), {
                         h_IM_4g, h_IM_gg, h_AntiPi0FitProb, h_AntiEtaFitProb, h_TreeFitProb
                     });
        return hists;
    }

    static cuttree::Cuts_t<Fill_t> GetCuts() {
        using cuttree::MultiCut_t;
        auto cuts = cuttree::ConvertCuts<Fill_t, CommonHist_t::Fill_t>(CommonHist_t::GetCuts());

        // Goldhaber cuts reduce pi0pi0 and pi0eta backgrounds
        auto goldhaber_cut = [] (const Fill_t& f) {
            const auto& tree = f.Shared;
            const double pi0 = ParticleTypeDatabase::Pi0.Mass();
            const double eta = ParticleTypeDatabase::Eta.Mass();
            const vec2 Pi0Pi0(pi0, pi0);
            const vec2 EtaPi0(eta, pi0);

            auto check_within = [] (vec2 im, vec2 center, double radius, double scale) {
                vec2 diff(im-center);
                diff.y *= scale;
                return diff.R() < radius;
            };

            for(unsigned i=0;i<tree.gg_gg1().size();i++) {
                const double im1 = tree.gg_gg1()[i];
                const double im2 = tree.gg_gg2()[i];
                if(   im1 < pi0+20
                   && im2 < pi0+20)
                    return false;
                if(check_within({im1, im2}, Pi0Pi0, 40, 1.0))
                    return false;
                if(check_within({im1, im2}, EtaPi0, 40, 1.5))
                    return false;
                if(check_within({im2, im1}, EtaPi0, 40, 1.5))
                    return false;
            }
            return true;
        };

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"AntiPi0FitProb<0.002", [] (const Fill_t& f) { return f.Shared.AntiPi0FitProb<0.002; } },
                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"AntiEtaFitProb<0.005", [] (const Fill_t& f) { return f.Shared.AntiEtaFitProb<0.005; } },
                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"Goldhaber", goldhaber_cut },
                              {"IM_gg", [] (const Fill_t& f) { return 180<f.Tree.IM_gg && f.Tree.IM_gg<460; } },
                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"TreeFitProb>0.04", [] (const Fill_t& f) { return f.Tree.TreeFitProb>0.04; } },
                              {"TreeFitProb>0.01", [] (const Fill_t& f) { return f.Tree.TreeFitProb>0.01; } },
                          });


        return cuts;
    }
};

struct SigPi0Hist_t : SigHist_t {
    using Tree_t = physics::EtapOmegaG::Sig_t::Pi0_t::BaseTree_t;

    TH2D* h_IM_3g_4g_low;     // Omega IM vs. EtaPrime IM
    TH2D* h_IM_3g_4g_high;     // Omega IM vs. EtaPrime IM

    struct Fill_t : SigHist_t::Fill_t {
        const Tree_t& Pi0;
        Fill_t(const CommonHist_t::Tree_t& common,
               const SharedTree_t& shared,
               const Tree_t& pi0) :
            SigHist_t::Fill_t(common, shared, pi0),
            Pi0(pi0)
        {}
    };

    SigPi0Hist_t(HistogramFactory HistFac) : SigHist_t(HistFac) {
        h_IM_3g_4g_low = HistFac.makeTH2D("Best #omega vs. #eta' IM",
                                          "IM(#pi^{0}#gamma#gamma) / MeV",
                                          "IM(#pi^{0}#gamma) / MeV",
                                          bins_IM_Etap, bins_IM_Omega,"h_IM_3g_4g_low"
                                          );
        h_IM_3g_4g_high = HistFac.makeTH2D("Best #omega vs. #eta' IM",
                                           "IM(#pi^{0}#gamma#gamma) / MeV",
                                           "IM(#pi^{0}#gamma) / MeV",
                                           bins_IM_Etap, bins_IM_Omega,"h_IM_3g_4g_high"
                                           );
    }

    void Fill(const Fill_t& f) const {
        SigHist_t::Fill(f);
        const Tree_t& pi0 = f.Pi0;
        h_IM_3g_4g_low->Fill(pi0.IM_Pi0gg, pi0.IM_Pi0g()[0], f.TaggW());
        h_IM_3g_4g_high->Fill(pi0.IM_Pi0gg, pi0.IM_Pi0g()[1], f.TaggW());
    }

    static cuttree::Cuts_t<Fill_t> GetCuts() {
        using cuttree::MultiCut_t;
        auto cuts = cuttree::ConvertCuts<Fill_t, SigHist_t::Fill_t>(SigHist_t::GetCuts());
        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"IM_Pi0g[1]>650", [] (const Fill_t& f) {
                                   return f.Pi0.IM_Pi0g()[1] > 650;
                               }},
                          });
        return cuts;
    }

};

struct SigOmegaPi0Hist_t : SigHist_t {
    using Tree_t = physics::EtapOmegaG::Sig_t::OmegaPi0_t::BaseTree_t;
    struct Fill_t : SigHist_t::Fill_t {
        const Tree_t& OmegaPi0;
        Fill_t(const CommonHist_t::Tree_t& common,
               const SharedTree_t& shared,
               const Tree_t& omegapi0) :
            SigHist_t::Fill_t(common, shared, omegapi0),
            OmegaPi0(omegapi0)
        {}
    };

    TH1D* h_Bachelor_E;

    SigOmegaPi0Hist_t(HistogramFactory HistFac) : SigHist_t(HistFac) {
        h_Bachelor_E = HistFac.makeTH1D("E_#gamma in #eta' frame","E_{#gamma} / MeV","",
                                        BinSettings(150,0,350),"h_Bachelor_E");

    }

    std::vector<TH1*> GetHists() const {
        auto hists = SigHist_t::GetHists();
        hists.insert(hists.end(), {h_Bachelor_E});
        return hists;
    }

    void Fill(const Fill_t& f) const {
        SigHist_t::Fill(f);
        const Tree_t& omegapi0 = f.OmegaPi0;
        h_Bachelor_E->Fill(omegapi0.Bachelor_E, f.TaggW());
    }

    static cuttree::Cuts_t<Fill_t> GetCuts() {
        return cuttree::ConvertCuts<Fill_t, SigHist_t::Fill_t>(SigHist_t::GetCuts());
    }

};

struct RefHist_t : CommonHist_t {
    using Tree_t = physics::EtapOmegaG::Ref_t::Tree_t;

    struct Fill_t : CommonHist_t::Fill_t {
        const Tree_t& Tree;
        Fill_t(const CommonHist_t::Tree_t& common, const Tree_t& tree) :
            CommonHist_t::Fill_t(common),
            Tree(tree)
        {}
    };

    TH1D* h_IM_2g;

    RefHist_t(HistogramFactory HistFac) : CommonHist_t(HistFac) {
        h_IM_2g = HistFac.makeTH1D("IM 2g","IM / MeV","",BinSettings(1100,0,1100),"h_IM_2g");
    }

    void Fill(const Fill_t& f) const {
        CommonHist_t::Fill(f);
        const Tree_t& tree = f.Tree;
        h_IM_2g->Fill(tree.IM_2g, f.TaggW());
    }

    std::vector<TH1*> GetHists() const {
        auto hists = CommonHist_t::GetHists();
        hists.insert(hists.end(), {h_IM_2g});
        return hists;
    }

    static cuttree::Cuts_t<Fill_t> GetCuts() {
        return cuttree::ConvertCuts<Fill_t, CommonHist_t::Fill_t>(CommonHist_t::GetCuts());
    }
};

template<typename Hist_t>
cuttree::Tree_t<MCTrue_Splitter<Hist_t>> makeMCSplitTree(const HistogramFactory& HistFac,
                                                         const std::string& treename)
{
    return cuttree::Make<MCTrue_Splitter<Hist_t>>(HistFac,
                                                  treename,
                                                  Hist_t::GetCuts()
                                                  );
}

int main(int argc, char** argv) {
    SetupLogger();

    signal(SIGINT, [] (int) { interrupt = true; } );

    TCLAP::CmdLine cmd("plot", ' ', "0.1");
    auto cmd_input = cmd.add<TCLAP::ValueArg<string>>("i","input","Input file",true,"","input");
    auto cmd_batchmode = cmd.add<TCLAP::MultiSwitchArg>("b","batch","Run in batch mode (no ROOT shell afterwards)",false);
    auto cmd_maxevents = cmd.add<TCLAP::MultiArg<int>>("m","maxevents","Process only max events",false,"maxevents");
    auto cmd_output = cmd.add<TCLAP::ValueArg<string>>("o","output","Output file",false,"","filename");

    auto cmd_tree = cmd.add<TCLAP::ValueArg<string>>("t","tree","Tree name",false,"Fitted/SigAll","treename");

    cmd.parse(argc, argv);

    WrapTFileInput input(cmd_input->getValue());

    auto link_branches = [&input] (const string treename, WrapTTree* wraptree, long long expected_entries) {
        TTree* t;
        if(!input.GetObject(treename,t))
            throw runtime_error("Cannot find tree "+treename+" in input file");
        if(expected_entries>=0 && t->GetEntries() != expected_entries)
            throw runtime_error("Tree "+treename+" does not have entries=="+to_string(expected_entries));
        if(wraptree->Matches(t)) {
            wraptree->LinkBranches(t);
            return true;
        }
        return false;
    };


    CommonHist_t::Tree_t treeCommon;
    if(!link_branches("EtapOmegaG/treeCommon", addressof(treeCommon), -1)) {
        LOG(ERROR) << "Cannot link branches of treeCommon";
        return 1;
    }

    auto entries = treeCommon.Tree->GetEntries();



    SigPi0Hist_t::Tree_t treeSigPi0;
    SigOmegaPi0Hist_t::Tree_t treeSigOmegaPi0;
    RefHist_t::Tree_t treeRef;

    // auto-detect which tree "type" to use
    const auto& treename = cmd_tree->getValue();
    if(link_branches("EtapOmegaG/"+treename, addressof(treeSigPi0), entries)) {
        LOG(INFO) << "Identified " << treename << " as signal tree (Pi0)";
    }
    else if(link_branches("EtapOmegaG/"+treename, addressof(treeSigOmegaPi0), entries)) {
        LOG(INFO) << "Identified " << treename << " as signal tree (OmegaPi0)";
    }
    else if(link_branches("EtapOmegaG/"+treename, addressof(treeRef), entries)) {
        LOG(INFO) << "Identified " << treename << " as reference tree";
    }
    else {
        LOG(ERROR) << "Could not identify " << treename;
        return 1;
    }

    SigHist_t::SharedTree_t treeSigShared;
    if(treeSigPi0 || treeSigOmegaPi0) {
        if(!link_branches("EtapOmegaG/SigShared", addressof(treeSigShared), entries)) {
            LOG(ERROR) << "Cannot find SigShared tree";
            return 1;
        }
    }

    unique_ptr<WrapTFileOutput> masterFile;
    if(cmd_output->isSet()) {
        masterFile = std_ext::make_unique<WrapTFileOutput>(cmd_output->getValue(),
                                                    WrapTFileOutput::mode_t::recreate,
                                                     true); // cd into masterFile upon creation
    }

    using MCSigPi0Hist_t = MCTrue_Splitter<SigPi0Hist_t>;
    using MCSigOmegaPi0Hist_t = MCTrue_Splitter<SigOmegaPi0Hist_t>;
    using MCRefHist_t = MCTrue_Splitter<RefHist_t>;

    HistogramFactory HistFac("EtapOmegaG");

    const auto& sanitized_treename = std_ext::replace_str(cmd_tree->getValue(),"/","_");
    auto cuttreeSigPi0 = treeSigPi0 ? makeMCSplitTree<SigPi0Hist_t>(HistFac, sanitized_treename) : nullptr;
    auto cuttreeSigOmegaPi0 = treeSigOmegaPi0 ? makeMCSplitTree<SigOmegaPi0Hist_t>(HistFac, sanitized_treename) : nullptr;
    auto cuttreeRef = treeRef ? makeMCSplitTree<RefHist_t>(HistFac, sanitized_treename) : nullptr;

    LOG(INFO) << "Tree entries=" << entries;
    auto max_entries = entries;
    if(cmd_maxevents->isSet() && cmd_maxevents->getValue().back()<entries) {
        max_entries = cmd_maxevents->getValue().back();
        LOG(INFO) << "Running until " << max_entries;
    }

    for(long long entry=0;entry<max_entries;entry++) {
        if(interrupt)
            break;

        treeCommon.Tree->GetEntry(entry);

        // we handle the Ref/Sig cut here to save some reading work
        if(treeSigShared && treeCommon.IsSignal) {
            treeSigShared.Tree->GetEntry(entry);
            if(treeSigPi0) {
                treeSigPi0.Tree->GetEntry(entry);
                cuttree::Fill<MCSigPi0Hist_t>(cuttreeSigPi0, {treeCommon, treeSigShared, treeSigPi0});
            }
            else if(treeSigOmegaPi0) {
                treeSigOmegaPi0.Tree->GetEntry(entry);
                cuttree::Fill<MCSigOmegaPi0Hist_t>(cuttreeSigOmegaPi0, {treeCommon, treeSigShared, treeSigOmegaPi0});
            }
        }
        else if(treeRef && !treeCommon.IsSignal) {
            treeRef.Tree->GetEntry(entry);
            cuttree::Fill<MCRefHist_t>(cuttreeRef, {treeCommon, treeRef});
        }
        if(entry % 100000 == 0)
            LOG(INFO) << "Processed " << 100.0*entry/entries << " %";
    }

    if(!cmd_batchmode->isSet()) {
        if(!std_ext::system::isInteractive()) {
            LOG(INFO) << "No TTY attached. Not starting ROOT shell.";
        }
        else {

            TRint app("EtapOmegaG_plot",&argc,argv,nullptr,0,true);

            if(masterFile)
                LOG(INFO) << "Stopped running, but close ROOT properly to write data to disk.";

            app.Run(kTRUE); // really important to return...
            if(masterFile)
                LOG(INFO) << "Writing output file...";
            masterFile = nullptr;   // and to destroy the master WrapTFile before TRint is destroyed
        }
    }


    return 0;
}