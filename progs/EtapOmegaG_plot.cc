#include "base/Logger.h"

#include "TRint.h"
#include "analysis/plot/root_draw.h"
#include "analysis/physics/etaprime/etaprime_omega_gamma.h"

#include "base/CmdLine.h"
#include "base/interval.h"

#include "base/printable.h"
#include "base/iterators.h"
#include "base/std_ext/string.h"
#include "base/WrapTFile.h"

#include "TH1D.h"

#include <functional>
#include <string>
#include <list>

using namespace ant;
using namespace ant::analysis;
using namespace std;


using CommonTree_t = physics::EtapOmegaG::TreeCommon;

template<typename SigRefTree_t>
struct Cut_t {
    using Passes_t = function<bool(const CommonTree_t&, const SigRefTree_t&)>;
    Cut_t(const std::string& name,
          Passes_t passes = [] (const CommonTree_t&, const SigRefTree_t&) { return true; }) : Name(name), Passes(passes) {}
    const std::string Name;
    const Passes_t Passes;
};


template<typename SigRefTree_t, typename Hist_t>
struct Node_t {
    SmartHistFactory HistFac; // directory for cuts
    SmartHistFactory H;       // subdir for hists only
    map<unsigned, Hist_t> Hists;
    typename Cut_t<SigRefTree_t>::Passes_t PassesCut;

    static map<unsigned, Hist_t> MakeHists(SmartHistFactory& h, const std::string& prefix)
    {
        h.SetTitlePrefix(prefix);
        map<unsigned, Hist_t> m;
        // 0=Data, 1=Sig, 2=Ref, 3=AllMC, 9=Unknown_Bkg, >=10 Known Bkg
        // create them here already to get them in the right order
        m.emplace(0, SmartHistFactory("Data", h, "Data"));
        m.emplace(3, SmartHistFactory("MC",   h, "MC"));
        m.emplace(1, SmartHistFactory("Sig",  h, "Sig"));
        m.emplace(2, SmartHistFactory("Ref",  h, "Ref"));
        return m;
    }

    Node_t(const SmartHistFactory& parentHistFac, Cut_t<SigRefTree_t> cut) :
        HistFac(cut.Name, parentHistFac, cut.Name),
        H("h", HistFac),
        Hists(MakeHists(H, HistFac.GetTitlePrefix())),
        PassesCut(cut.Passes)
    {}

    bool Fill(const CommonTree_t& treeCommon, const SigRefTree_t& treeSigRef) {
       if(PassesCut(treeCommon, treeSigRef)) {

           const auto mctrue = treeCommon.MCTrue();
           auto it_hist = Hists.lower_bound(mctrue);
           if(it_hist == Hists.end() || it_hist->first != mctrue) {
               // mctrue should be 9 or higher here...
               const string& name = mctrue>=10 ?
                                        physics::EtapOmegaG::ptreeBackgrounds[mctrue-10].Name
                                    : "Other";
               it_hist = Hists.emplace_hint(it_hist, mctrue, SmartHistFactory("Bkg_"+name, H, name));
           }

           it_hist->second.Fill(treeCommon, treeSigRef);
           if(mctrue>0) // mctrue is never 3 in tree, use this to sum up all MC
               Hists.at(3).Fill(treeCommon, treeSigRef);
           return true;
       }
       return false;
    }
};


template<typename SigRefTree_t>
using MultiCut_t = std::vector<Cut_t<SigRefTree_t>>;

template<typename SigRefTree_t>
using Cuts_t = std::list<MultiCut_t<SigRefTree_t>>;

template<typename SigRefTree_t>
using CutsIterator_t = typename Cuts_t<SigRefTree_t>::const_iterator;


template<typename SigRefTree_t, typename Hist_t>
using CutTree_t = typename Tree<Node_t<SigRefTree_t, Hist_t>>::node_t;

template<typename SigRefTree_t, typename Hist_t>
void BuildCutTree(CutTree_t<SigRefTree_t, Hist_t> cuttree, CutsIterator_t<SigRefTree_t> first, CutsIterator_t<SigRefTree_t> last) {
    if(first == last)
        return;
    const MultiCut_t<SigRefTree_t>& multicut = *first;
    for(const Cut_t<SigRefTree_t>& cut : multicut) {
        auto daughter = cuttree->CreateDaughter(cuttree->Get().HistFac, cut);
        BuildCutTree<SigRefTree_t, Hist_t>(daughter, std::next(first), last);
    }
}

template<typename Hist_t, typename SigRefTree_t = typename Hist_t::Tree_t>
CutTree_t<SigRefTree_t, Hist_t> MakeCutTree(SmartHistFactory histFac, const string& name) {
    const auto& cuts = Hist_t::GetCuts();
    auto cuttree = Tree<Node_t<SigRefTree_t, Hist_t>>::MakeNode(histFac, Cut_t<SigRefTree_t>{name});
    BuildCutTree<SigRefTree_t, Hist_t>(cuttree, cuts.begin(), cuts.end());
    return cuttree;
}

template<typename Hist_t, typename SigRefTree_t = typename Hist_t::Tree_t>
void FillCutTree(CutTree_t<SigRefTree_t, Hist_t> cuttree, const CommonTree_t& commonTree, const SigRefTree_t& tree) {
    if(cuttree->Get().Fill(commonTree, tree)) {
        for(const auto& d : cuttree->Daughters()) {
            FillCutTree<Hist_t, SigRefTree_t>(d, commonTree, tree);
        }
    }
}

// define the structs containing the histograms
// and the cuts. for simple branch variables, that could
// be combined...


struct CommonHist_t {
    TH1D* h_KinFitChi2;

    CommonHist_t(SmartHistFactory HistFac) {
        h_KinFitChi2 = HistFac.makeTH1D("KinFitChi2","#chi^{2}","",BinSettings(200,0,100),"h_KinFitChi2");
    }
    void Fill(const CommonTree_t& treeCommon) {
        h_KinFitChi2->Fill(treeCommon.KinFitChi2, treeCommon.TaggW);
    }

    // there are some common cuts for Sig and Ref channel
    template<typename Tree_t>
    static Cuts_t<Tree_t> GetCuts() {
        Cuts_t<Tree_t> cuts;
        cuts.emplace_back(MultiCut_t<Tree_t>{
                                 {"CBSumVeto=0", [] (const CommonTree_t& treeCommon, const Tree_t&) { return treeCommon.CBSumVetoE==0; } },
                                 {"CBSumVeto<0.25", [] (const CommonTree_t& treeCommon, const Tree_t&) { return treeCommon.CBSumVetoE<0.25; } },
                                 {"PIDSumE=0", [] (const CommonTree_t& treeCommon, const Tree_t&) { return treeCommon.PIDSumE==0; } },
                                 {"PIDSumE<0.25", [] (const CommonTree_t& treeCommon, const Tree_t&) { return treeCommon.PIDSumE<0.25; } },
                             });
        cuts.emplace_back(MultiCut_t<Tree_t>{
                                 {"KinFitChi2<10", [] (const CommonTree_t& treeCommon, const Tree_t&) { return treeCommon.KinFitChi2<10; } },
                                 {"KinFitChi2<20", [] (const CommonTree_t& treeCommon, const Tree_t&) { return treeCommon.KinFitChi2<20; } },
                             });
        return cuts;
    }

};

struct SigHist_t : CommonHist_t {
    using Tree_t = physics::EtapOmegaG::Sig_t::Tree_t;

    TH1D* h_TreeFitChi2;
    TH1D* h_Bachelor_E;

    SigHist_t(SmartHistFactory HistFac) : CommonHist_t(HistFac) {
        h_TreeFitChi2 = HistFac.makeTH1D("TreeFitChi2", "#chi^{2}","",BinSettings(200,0,100),"h_TreeFitChi2");
        h_Bachelor_E = HistFac.makeTH1D("E_#gamma in #eta' frame","E_{#gamma}","",BinSettings(400,0,400),"h_Bachelor_E");
    }

    void Fill(const CommonTree_t& treeCommon, const Tree_t& treeSig) {
        CommonHist_t::Fill(treeCommon);
        h_TreeFitChi2->Fill(treeSig.TreeFitChi2, treeCommon.TaggW);
        h_Bachelor_E->Fill(treeSig.Bachelor_best_best, treeCommon.TaggW);
    }

    static Cuts_t<Tree_t> GetCuts() {
        auto cuts = CommonHist_t::GetCuts<Tree_t>();
        cuts.emplace_back(MultiCut_t<Tree_t>{
                              {"TreeFitChi2<20", [] (const CommonTree_t&, const Tree_t& tree) { return tree.TreeFitChi2<20; } },
                              {"TreeFitChi2<50", [] (const CommonTree_t&, const Tree_t& tree) { return tree.TreeFitChi2<50; } },
                          });
        return cuts;
    }
};

struct RefHist_t : CommonHist_t {
    using Tree_t = physics::EtapOmegaG::Ref_t::Tree_t;

    TH1D* h_IM_2g;

    RefHist_t(SmartHistFactory HistFac) : CommonHist_t(HistFac) {
        h_IM_2g = HistFac.makeTH1D("IM 2g","IM / MeV","",BinSettings(1100,0,1100),"h_IM_2g");
    }

    void Fill(const CommonTree_t& treeCommon, const Tree_t& treeRef) {
        CommonHist_t::Fill(treeCommon);
        h_IM_2g->Fill(treeRef.IM_2g, treeCommon.TaggW);
    }

    static Cuts_t<Tree_t> GetCuts() {
        // reference does not have more to cut than common stuff...
        return CommonHist_t::GetCuts<Tree_t>();
    }
};

int main(int argc, char** argv) {
    SetupLogger();

    TCLAP::CmdLine cmd("plot", ' ', "0.1");
    auto cmd_file = cmd.add<TCLAP::ValueArg<string>>("i","input","Input file",true,"","input");;
    cmd.parse(argc, argv);

    WrapTFileInput input(cmd_file->getValue());

    auto link_branches = [&input] (const string treename, WrapTTree* wraptree, long long expected_entries) {
        TTree* t;
        if(!input.GetObject(treename,t))
            throw runtime_error("Cannot find tree "+treename+" in input file");
        if(expected_entries>=0 && t->GetEntries() != expected_entries)
            throw runtime_error("Tree "+treename+" does not have entries=="+to_string(expected_entries));
        wraptree->LinkBranches(t);
    };


    CommonTree_t treeCommon;
    link_branches("EtapOmegaG/treeCommon", addressof(treeCommon), -1);
    const auto entries = treeCommon.Tree->GetEntries();

    SigHist_t::Tree_t treeSig;
    link_branches("EtapOmegaG/treeSig", addressof(treeSig), entries);
    SigHist_t::Tree_t treeSigFitted;
    link_branches("EtapOmegaG/treeSigFitted", addressof(treeSigFitted), entries);
    RefHist_t::Tree_t treeRef;
    link_branches("EtapOmegaG/treeRef", addressof(treeRef), entries);
    RefHist_t::Tree_t treeRefFitted;
    link_branches("EtapOmegaG/treeRefFitted", addressof(treeRefFitted), entries);

    LOG(INFO) << "Tree entries=" << entries;

    SmartHistFactory HistFac("EtapOmegaG");


    auto cuttreeSig = MakeCutTree<SigHist_t>(HistFac, "Sig");
    auto cuttreeSigFitted = MakeCutTree<SigHist_t>(HistFac, "SigFitted");
    auto cuttreeRef = MakeCutTree<RefHist_t>(HistFac, "Ref");
    auto cuttreeRefFitted = MakeCutTree<RefHist_t>(HistFac, "RefFitted");

    for(long long entry=0;entry<entries;entry++) {
        treeCommon.Tree->GetEntry(entry);

        // we handle the Ref/Sig cut here to save some reading work
        if(treeCommon.IsSignal) {
            treeSig.Tree->GetEntry(entry);
            FillCutTree<SigHist_t>(cuttreeSig, treeCommon, treeSig);
            treeSigFitted.Tree->GetEntry(entry);
            FillCutTree<SigHist_t>(cuttreeSigFitted, treeCommon, treeSigFitted);
        }
        else {
            treeRef.Tree->GetEntry(entry);
            FillCutTree<RefHist_t>(cuttreeRef, treeCommon, treeRef);
            treeRefFitted.Tree->GetEntry(entry);
            FillCutTree<RefHist_t>(cuttreeRefFitted, treeCommon, treeRefFitted);
        }
    }


    int fake_argc=0;
    char** fake_argv=nullptr;
    TRint app("plot",&fake_argc,fake_argv,nullptr,0,true);

    //canvas() << padoption::LogY << h_IM_2g << endc;

    app.Run();



    return 0;
}