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

class MyTInterruptHandler : public TSignalHandler {
public:
    MyTInterruptHandler() : TSignalHandler(kSigInterrupt, kFALSE) { }

    Bool_t  Notify() {
        if (fDelay) {
            fDelay++;
            return kTRUE;
        }
        interrupt = true;
        cout << " >>> Interrupted! " << endl;
        return kTRUE;
    }
};

template<typename Hist_t>
struct MCTrue_Splitter : cuttree::StackedHists_t<Hist_t> {

    // Hist_t should have that type defined
    using Fill_t = typename Hist_t::Fill_t;

    MCTrue_Splitter(const HistogramFactory& histFac) : cuttree::StackedHists_t<Hist_t>(histFac) {
        using cuttree::HistMod_t;
        this->GetHist(0, "Data", HistMod_t::MakeDataPoints(kBlack));
        this->GetHist(1, "Sig",  HistMod_t::MakeLine(kRed, 2.0));
        this->GetHist(2, "Ref",  HistMod_t::MakeLine(kRed, 2.0));
        // mctrue is never >=3 (and <9) in tree, use this to sum up all MC and all bkg MC
        // see also Fill()
        this->GetHist(3, "Sum_MC", HistMod_t::MakeLine(kBlack, 2.0));
        this->GetHist(4, "Bkg_MC", HistMod_t::MakeLine(kGray, 2.0));
    }

    void Fill(const Fill_t& f) {

        const unsigned mctrue = f.Common.MCTrue;

        auto get_bkg_name = [] (unsigned mctrue) {
            const string& name = mctrue>=10 ?
                                     physics::EtapOmegaG::ptreeBackgrounds[mctrue-10].Name
                                 : "Other";
            return "Bkg_"+name;
        };

        using cuttree::HistMod_t;
        const Hist_t& hist = mctrue<9 ? this->GetHist(mctrue) :
                                        this->GetHist(mctrue,
                                                      get_bkg_name(mctrue),
                                                      HistMod_t::MakeLine(HistMod_t::GetColor(mctrue-9))
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

    template<typename SigRefTree_t>
    struct SigRefFill_t : Fill_t {
        const SigRefTree_t& Tree;
        SigRefFill_t(const CommonHist_t::Tree_t& common, const SigRefTree_t& tree) :
            CommonHist_t::Fill_t(common),
            Tree(tree)
        {}
    };

    CommonHist_t(HistogramFactory) {}
    void Fill(const Fill_t&) const {}
    std::vector<TH1*> GetHists() const { return {}; }

    // Sig and Ref channel share some cuts...
    static cuttree::Cuts_t<Fill_t> GetCuts() {
        using cuttree::MultiCut_t;
        cuttree::Cuts_t<Fill_t> cuts;
        cuts.emplace_back(MultiCut_t<Fill_t>{
                                 {"CBSumVeto=0", [] (const Fill_t& f) { return f.Common.CBSumVetoE==0; } },
                                 {"CBSumVeto<0.25", [] (const Fill_t& f) { return f.Common.CBSumVetoE<0.25; } },
                                 {"PIDSumE=0", [] (const Fill_t& f) { return f.Common.PIDSumE==0; } },
                                 {"PIDSumE<0.25", [] (const Fill_t& f) { return f.Common.PIDSumE<0.25; } },
                             });
        return cuts;
    }


};

struct SigHist_t : CommonHist_t {
    using Tree_t = physics::EtapOmegaG::Sig_t::Tree_t;
    using Fill_t = CommonHist_t::SigRefFill_t<Tree_t>;

    TH2D* h_IM_gg_gg;     // Goldhaber plot
    TH1D* h_IM_4g;        // EtaPrime IM
    TH2D* h_IM_3g_4g;     // Omega IM vs. EtaPrime IM
    TH1D* h_TreeFitChi2;
    TH1D* h_Bachelor_E;

    SigHist_t(HistogramFactory HistFac) : CommonHist_t(HistFac) {
        BinSettings bins_goldhaber(400, 0, 900);
        const string axislabel_goldhaber("2#gamma IM / MeV");

        h_IM_gg_gg = HistFac.makeTH2D("IM 2#gamma-2#gamma",
                                    axislabel_goldhaber, axislabel_goldhaber,
                                    bins_goldhaber, bins_goldhaber,
                                    "h_IM_gg_gg"
                                    );
        BinSettings IM_Etap(100,800,1000);
        BinSettings IM_Omega(100, 700, 950);
        h_IM_4g = HistFac.makeTH1D("#eta' IM", "IM(#pi^{0}#gamma#gamma) / MeV","",IM_Etap,"h_IM_4g");
        h_IM_3g_4g = HistFac.makeTH2D("Best #omega vs. #eta' IM",
                                      "IM(#pi^{0}#gamma#gamma) / MeV",
                                      "IM(#pi^{0}#gamma) best / MeV",
                                      IM_Etap, IM_Omega,"h_IM_3g_4g"
                                      );
        h_TreeFitChi2 = HistFac.makeTH1D("TreeFitChi2", "#chi^{2}","",BinSettings(200,0,100),"h_TreeFitChi2");
        h_Bachelor_E = HistFac.makeTH1D("E_#gamma in #eta' frame","E_{#gamma} / MeV","",BinSettings(400,0,400),"h_Bachelor_E");
    }

    void Fill(const Fill_t& f) const {
        CommonHist_t::Fill(f);
        const Tree_t& tree = f.Tree;
        for(unsigned i=0;i<tree.gg_gg1().size();i++){
            h_IM_gg_gg->Fill(tree.gg_gg1()[i], tree.gg_gg2()[i], f.TaggW());
            h_IM_gg_gg->Fill(tree.gg_gg2()[i], tree.gg_gg1()[i], f.TaggW());
        }
        h_TreeFitChi2->Fill(tree.TreeFitChi2, f.TaggW());
        h_Bachelor_E->Fill(tree.Bachelor_E, f.TaggW());
        h_IM_3g_4g->Fill(tree.IM_Pi0gg_best, tree.IM_Pi0g_best, f.TaggW());
        h_IM_4g->Fill(tree.IM_Pi0gg_best, f.TaggW());
    }

    std::vector<TH1*> GetHists() const {
        auto hists = CommonHist_t::GetHists();
        hists.insert(hists.end(), {h_IM_4g, h_TreeFitChi2, h_Bachelor_E});
        return hists;
    }

    static cuttree::Cuts_t<Fill_t> GetCuts() {
        using cuttree::MultiCut_t;
        auto cuts = cuttree::ConvertCuts<Fill_t, CommonHist_t::Fill_t>(CommonHist_t::GetCuts());

        // reduces pi0pi0 and pi0eta backgrounds
        auto goldhaber_cut = [] (const Fill_t& f) {
            const Tree_t& tree = f.Tree;
            const double pi0 = ParticleTypeDatabase::Pi0.Mass();
            const double eta = ParticleTypeDatabase::Eta.Mass();
            const TVector2 Pi0Pi0(pi0, pi0);
            const TVector2 EtaPi0(eta, pi0);

            auto check_within = [] (TVector2 im, TVector2 center, double radius, double scale) {
                TVector2 diff(im-center);
                diff.Set(diff.X(), scale*diff.Y());
                return diff.Mod() < radius;
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
                              {"Goldhaber", goldhaber_cut },
                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"TreeFitChi2<20", [] (const Fill_t& f) { return f.Tree.TreeFitChi2<20; } },
                              {"TreeFitChi2<50", [] (const Fill_t& f) { return f.Tree.TreeFitChi2<50; } },
                          });


        return cuts;
    }
};

struct RefHist_t : CommonHist_t {
    using Tree_t = physics::EtapOmegaG::Ref_t::Tree_t;
    using Fill_t = CommonHist_t::SigRefFill_t<Tree_t>;

    TH1D* h_KinFitChi2;
    TH1D* h_IM_2g;

    RefHist_t(HistogramFactory HistFac) : CommonHist_t(HistFac) {
        h_KinFitChi2 = HistFac.makeTH1D("KinFitChi2","#chi^{2}","",BinSettings(200,0,100),"h_KinFitChi2");
        h_IM_2g = HistFac.makeTH1D("IM 2g","IM / MeV","",BinSettings(1100,0,1100),"h_IM_2g");
    }

    void Fill(const Fill_t& f) const {
        CommonHist_t::Fill(f);
        const Tree_t& tree = f.Tree;
        h_KinFitChi2->Fill(tree.KinFitChi2, f.TaggW());
        h_IM_2g->Fill(tree.IM_2g, f.TaggW());
    }

    std::vector<TH1*> GetHists() const {
        auto hists = CommonHist_t::GetHists();
        hists.insert(hists.end(), {h_KinFitChi2, h_IM_2g});
        return hists;
    }

    static cuttree::Cuts_t<Fill_t> GetCuts() {
        using cuttree::MultiCut_t;
        auto cuts = cuttree::ConvertCuts<Fill_t, CommonHist_t::Fill_t>(CommonHist_t::GetCuts());
        cuts.emplace_back(MultiCut_t<Fill_t>{
                                 {"KinFitChi2<10", [] (const Fill_t& f) { return f.Tree.KinFitChi2<10; } },
                                 {"KinFitChi2<20", [] (const Fill_t& f) { return f.Tree.KinFitChi2<20; } },
                             });
        return cuts;
    }
};

int main(int argc, char** argv) {
    SetupLogger();

    TCLAP::CmdLine cmd("plot", ' ', "0.1");
    auto cmd_input = cmd.add<TCLAP::ValueArg<string>>("i","input","Input file",true,"","input");
    auto cmd_batchmode = cmd.add<TCLAP::MultiSwitchArg>("b","batch","Run in batch mode (no ROOT shell afterwards)",false);
    auto cmd_maxevents = cmd.add<TCLAP::MultiArg<int>>("m","maxevents","Process only max events",false,"maxevents");
    auto cmd_output = cmd.add<TCLAP::ValueArg<string>>("o","output","Output file",false,"","filename");

    auto cmd_tree = cmd.add<TCLAP::ValueArg<string>>("t","tree","Tree name",false,"Fitted/SigAll","treename");

    cmd.parse(argc, argv);

    int fake_argc=1;
    char* fake_argv[2];
    fake_argv[0] = argv[0];
    if(cmd_batchmode->isSet()) {
        fake_argv[fake_argc++] = strdup("-b");
    }
    TRint app("EtapOmegaG_plot",&fake_argc,fake_argv,nullptr,0,true);
    auto oldsig = app.GetSignalHandler();
    oldsig->Remove();
    auto mysig = new MyTInterruptHandler();
    mysig->Add();
    gSystem->AddSignalHandler(mysig);

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



    SigHist_t::Tree_t treeSig;
    RefHist_t::Tree_t treeRef;

    // auto-detect which tree "type" to use
    const auto& treename = cmd_tree->getValue();
    if(link_branches("EtapOmegaG/"+treename, addressof(treeSig), entries)) {
        LOG(INFO) << "Identified " << treename << " as signal tree";
    }
    else if(link_branches("EtapOmegaG/"+treename, addressof(treeRef), entries)) {
        LOG(INFO) << "Identified " << treename << " as reference tree";
    }
    else {
        LOG(ERROR) << "Could not identify " << treename;
        return 1;
    }

    unique_ptr<WrapTFileOutput> masterFile;
    if(cmd_output->isSet()) {
        masterFile = std_ext::make_unique<WrapTFileOutput>(cmd_output->getValue(),
                                                    WrapTFileOutput::mode_t::recreate,
                                                     true); // cd into masterFile upon creation
    }

    using MCSigHist_t = MCTrue_Splitter<SigHist_t>;
    using MCRefHist_t = MCTrue_Splitter<RefHist_t>;

    HistogramFactory HistFac("EtapOmegaG");

    const auto& sanitized_treename = std_ext::replace_str(cmd_tree->getValue(),"/","_");
    auto cuttreeSig = treeSig ? cuttree::Make<MCSigHist_t>(HistFac,
                                                           sanitized_treename,
                                                           SigHist_t::GetCuts()
                                                           ) : nullptr;
    auto cuttreeRef = treeRef ? cuttree::Make<MCRefHist_t>(HistFac,
                                                           sanitized_treename,
                                                           RefHist_t::GetCuts()
                                                           ) : nullptr;

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
        if(treeSig && treeCommon.IsSignal) {
            treeSig.Tree->GetEntry(entry);
            cuttree::Fill<MCSigHist_t>(cuttreeSig, {treeCommon, treeSig});
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

            mysig->Remove();
            oldsig->Add();
            gSystem->AddSignalHandler(oldsig);
            delete mysig;

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