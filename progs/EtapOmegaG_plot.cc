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
    SmartHistFactory HistFac;
    Hist_t Hist;
    typename Cut_t<SigRefTree_t>::Passes_t PassesCut;


    Node_t(const SmartHistFactory& parentHistFac, Cut_t<SigRefTree_t> cut) :
        HistFac(cut.Name, parentHistFac, cut.Name), Hist(HistFac), PassesCut(cut.Passes) {}

    bool Fill(const CommonTree_t& treeCommon, const SigRefTree_t& treeSigRef) {
       if(PassesCut(treeCommon, treeSigRef)) {
           Hist.Fill(treeCommon, treeSigRef);
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

template<typename SigRefTree_t, typename Hist_t>
void FillCutTree(CutTree_t<SigRefTree_t, Hist_t> cuttree, const CommonTree_t& commonTree, const SigRefTree_t& tree) {
    if(cuttree->Get().Fill(commonTree, tree)) {
        for(const auto& d : cuttree->Daughters()) {
            FillCutTree<SigRefTree_t, Hist_t>(d, commonTree, tree);
        }
    }
}

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

    using SigTree_t = physics::EtapOmegaG::Sig_t::Tree_t;
    using RefTree_t = physics::EtapOmegaG::Ref_t::Tree_t;

    SigTree_t treeSig;
    link_branches("EtapOmegaG/treeSig", addressof(treeSig), entries);
    SigTree_t treeSigFitted;
    link_branches("EtapOmegaG/treeSigFitted", addressof(treeSigFitted), entries);
    RefTree_t treeRef;
    link_branches("EtapOmegaG/treeRef", addressof(treeRef), entries);
    RefTree_t treeRefFitted;
    link_branches("EtapOmegaG/treeRefFitted", addressof(treeRefFitted), entries);

    LOG(INFO) << "Tree entries=" << entries;

    SmartHistFactory HistFac("EtapOmegaG");

    struct CommonHist_t {
        TH1D* h_KinFitChi2;
        CommonHist_t(SmartHistFactory HistFac) {
            h_KinFitChi2 = HistFac.makeTH1D("KinFitChi2","#chi^2","",BinSettings(200,0,100),"h_KinFitChi2");
        }
        void Fill(const CommonTree_t& treeCommon) {
            h_KinFitChi2->Fill(treeCommon.KinFitChi2, treeCommon.TaggW);
        }
    };

    struct RefHist_t : CommonHist_t {
        TH1D* h_IM_2g;
        RefHist_t(SmartHistFactory HistFac) : CommonHist_t(HistFac) {
            h_IM_2g = HistFac.makeTH1D("IM 2g","IM / MeV","",BinSettings(1100,0,1100),"h_IM_2g");
        }
        void Fill(const CommonTree_t& treeCommon, const RefTree_t& treeSigRef) {
            CommonHist_t::Fill(treeCommon);
            h_IM_2g->Fill(treeSigRef.IM_2g, treeCommon.TaggW);
        }
    };

    Cuts_t<RefTree_t> refCuts;
    refCuts.emplace_back(MultiCut_t<RefTree_t>{
                             {"CBSumVeto=0", [] (const CommonTree_t& treeCommon, const RefTree_t&) { return treeCommon.CBSumVetoE==0; } },
                             {"CBSumVeto<0.25", [] (const CommonTree_t& treeCommon, const RefTree_t&) { return treeCommon.CBSumVetoE<0.25; } },
                             {"PIDSumE=0", [] (const CommonTree_t& treeCommon, const RefTree_t&) { return treeCommon.PIDSumE==0; } },
                             {"PIDSumE<0.25", [] (const CommonTree_t& treeCommon, const RefTree_t&) { return treeCommon.PIDSumE<0.25; } },
                         });
    refCuts.emplace_back(MultiCut_t<RefTree_t>{
                             {"KinFitChi2<5", [] (const CommonTree_t& treeCommon, const RefTree_t&) { return treeCommon.KinFitChi2<10; } },
                             {"KinFitChi2<15", [] (const CommonTree_t& treeCommon, const RefTree_t&) { return treeCommon.KinFitChi2<100; } },
                         });

    auto cuttreeRef = Tree<Node_t<RefTree_t, RefHist_t>>::MakeNode(HistFac, Cut_t<RefTree_t>{"Ref"});
    BuildCutTree<RefTree_t, RefHist_t>(cuttreeRef, refCuts.begin(), refCuts.end());

    auto cuttreeRefFitted = Tree<Node_t<RefTree_t, RefHist_t>>::MakeNode(HistFac, Cut_t<RefTree_t>{"RefFitted"});
    BuildCutTree<RefTree_t, RefHist_t>(cuttreeRefFitted, refCuts.begin(), refCuts.end());

    for(long long entry=0;entry<entries;entry++) {
        treeCommon.Tree->GetEntry(entry);

        // we handle the Ref/Sig cut here to save some reading work
        if(treeCommon.IsSignal) {
            treeSig.Tree->GetEntry(entry);
            treeSigFitted.Tree->GetEntry(entry);
        }
        else {
            treeRef.Tree->GetEntry(entry);
            treeRefFitted.Tree->GetEntry(entry);
            FillCutTree<RefTree_t, RefHist_t>(cuttreeRef, treeCommon, treeRef);
            FillCutTree<RefTree_t, RefHist_t>(cuttreeRefFitted, treeCommon, treeRefFitted);
        }
    }


    int fake_argc=0;
    char** fake_argv=nullptr;
    TRint app("plot",&fake_argc,fake_argv,nullptr,0,true);

    //canvas() << padoption::LogY << h_IM_2g << endc;

    app.Run();



    return 0;
}