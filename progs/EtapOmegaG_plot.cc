#include "base/Logger.h"

#include "analysis/plot/CutTree.h"
#include "analysis/physics/etaprime/etaprime_omega_gamma.h"

#include "expconfig/ExpConfig.h"
#include "expconfig/detectors/EPT.h"

#include "base/CmdLine.h"
#include "base/interval.h"
#include "base/printable.h"
#include "base/WrapTFile.h"
#include "base/std_ext/string.h"
#include "base/std_ext/system.h"
#include "base/ProgressCounter.h"

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

    MCTrue_Splitter(const HistogramFactory& histFac,
                    const cuttree::TreeInfo_t& treeInfo) :
        cuttree::StackedHists_t<Hist_t>(histFac, treeInfo)
    {
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
    using Shared_t = physics::EtapOmegaG::SharedTree_t;
    struct Fill_t {
        const Tree_t& Common;
        const Shared_t& Shared;

        Fill_t(const Tree_t& common, const Shared_t& shared) :
            Common(common), Shared(shared) {}
        double TaggW() const {
            return Common.TaggW;
        }

        struct LogProb_t {
            LogProb_t(const double& prob_) : prob(prob_) {}
            operator double() const { return std::log10(prob); }
        private:
            const double& prob;
        };

        LogProb_t KinFitProb{Shared.KinFitProb};
    };



    const BinSettings bins_FitProb{100, -5, 0};
    TH1D* h_KinFitProb = nullptr;
    TH1D* h_CBSumE = nullptr;
    TH1D* h_CBSumVetoE = nullptr;
    TH1D* h_PIDSumE = nullptr;
    TH1D* h_MissingMass = nullptr;
    TH1D* h_DiscardedEk = nullptr;
    TH2D* h_ProtonTOF = nullptr;
    TH2D* h_ProtonTOFFitted = nullptr;
    TH2D* h_ProtonVetoE = nullptr;
    TH2D* h_ProtonShortE = nullptr;



    const bool isLeaf = false;


    CommonHist_t(HistogramFactory HistFac, cuttree::TreeInfo_t treeInfo) :
        isLeaf(treeInfo.nDaughters==0)
    {
        h_KinFitProb = HistFac.makeTH1D("KinFitProb","log p","",bins_FitProb,"h_KinFitProb");
        h_CBSumE = HistFac.makeTH1D("CB Sum E","E / MeV","",BinSettings(100,500,1600),"h_CBSumE");
        h_CBSumVetoE = HistFac.makeTH1D("CB Veto Sum E","E / MeV","",BinSettings(50,0,10),"h_CBSumVetoE");
        h_PIDSumE = HistFac.makeTH1D("PID Sum E","E / MeV","",BinSettings(50,0,10),"h_PIDSumE");
        h_MissingMass = HistFac.makeTH1D("MissingMass","m / MeV","",BinSettings(200,600,1300),"h_MissingMass");
        h_DiscardedEk = HistFac.makeTH1D("DiscardedEk","E / MeV","",BinSettings(200,0,500),"h_DiscardedEk");
        if(!isLeaf)
            return;
        BinSettings bins_protonE(100,0,600);
        h_ProtonTOF = HistFac.makeTH2D("ProtonTOF","t","E",
                                       BinSettings(50,-5,20),bins_protonE,"h_ProtonTOF");
        h_ProtonTOFFitted = HistFac.makeTH2D("ProtonTOFFitted","t","E fitted",
                                             BinSettings(50,-5,20),bins_protonE,"h_ProtonTOFFitted");
        h_ProtonVetoE = HistFac.makeTH2D("ProtonVeto","E fitted","Veto E",
                                        bins_protonE,BinSettings(50,0,8),"h_ProtonVeto");
        h_ProtonShortE = HistFac.makeTH2D("ProtonShortE","E fitted","E short",
                                          bins_protonE,BinSettings(100,0,300),"h_ProtonShortE");
    }


    void Fill(const Fill_t& f) const {

        h_KinFitProb->Fill(f.KinFitProb, f.TaggW());
        h_CBSumE->Fill(f.Common.CBSumE, f.TaggW());
        h_CBSumVetoE->Fill(f.Shared.CBSumVetoE, f.TaggW());
        h_PIDSumE->Fill(f.Common.PIDSumE, f.TaggW());
        h_MissingMass->Fill(f.Shared.MissingMass, f.TaggW());
        h_DiscardedEk->Fill(f.Shared.DiscardedEk, f.TaggW());
        if(!isLeaf)
            return;
        h_ProtonTOF->Fill(f.Shared.ProtonTime, f.Shared.ProtonE, f.TaggW());
        h_ProtonTOFFitted->Fill(f.Shared.ProtonTime, f.Shared.FittedProtonE, f.TaggW());
        h_ProtonVetoE->Fill(f.Shared.FittedProtonE, f.Shared.ProtonVetoE, f.TaggW());
        h_ProtonShortE->Fill(f.Shared.FittedProtonE, f.Shared.ProtonShortE, f.TaggW());
    }

    std::vector<TH1*> GetHists() const {
        return {h_KinFitProb, h_CBSumE, h_CBSumVetoE, h_PIDSumE, h_MissingMass,
                    h_DiscardedEk};
    }

    // Sig and Ref channel (can) share some cuts...
    static cuttree::Cuts_t<Fill_t> GetCuts() {
        using cuttree::MultiCut_t;
        cuttree::Cuts_t<Fill_t> cuts;
        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"DiscardedEk=0", [] (const Fill_t& f) { return f.Shared.DiscardedEk == 0; } },
                              {"DiscardedEk<50", [] (const Fill_t& f) { return f.Shared.DiscardedEk < 50; } },
                          });
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
            CommonHist_t::Fill_t(common, shared),
            Shared(shared),
            Tree(tree)
        {}

        LogProb_t AntiPi0FitProb{Shared.AntiPi0FitProb};
        LogProb_t AntiEtaFitProb{Shared.AntiEtaFitProb};
        LogProb_t TreeFitProb{Tree.TreeFitProb};
    };

    TH1D* h_IM_4g;        // EtaPrime IM

    TH1D* h_AntiPi0FitProb;
    TH1D* h_AntiEtaFitProb;
    TH1D* h_TreeFitProb;

    TH1D* h_AntiPi0ZVertex;
    TH1D* h_AntiEtaZVertex;
    TH1D* h_TreeZVertex;

    TH2D* h_Pi0EtaFitProb;
    TH2D* h_Pi0TreeFitProb;
    TH2D* h_EtaTreeFitProb;

    TH1D* h_Bachelor_E;

    const BinSettings bins_IM_Etap {100, 800,1050};
    const BinSettings bins_IM_Omega{100, 550, 950};
    const BinSettings bins_ZVertex{100, -15, 15};

    SigHist_t(HistogramFactory HistFac, cuttree::TreeInfo_t treeInfo) : CommonHist_t(HistFac, treeInfo) {
        h_IM_4g = HistFac.makeTH1D("#eta' IM", "IM(#pi^{0}#gamma#gamma) / MeV","",bins_IM_Etap,"h_IM_4g");

        h_AntiPi0FitProb = HistFac.makeTH1D("AntiPi0FitProb", "log p","",bins_FitProb,"h_AntiPi0FitProb");
        h_AntiEtaFitProb = HistFac.makeTH1D("AntiEtaFitProb", "log p","",bins_FitProb,"h_AntiEtaFitProb");
        h_TreeFitProb = HistFac.makeTH1D("TreeFitProb", "log p","",bins_FitProb,"h_TreeFitProb");

        h_AntiPi0ZVertex = HistFac.makeTH1D("AntiPi0ZVertex", "z / cm","",bins_ZVertex,"h_AntiPi0ZVertex");
        h_AntiEtaZVertex = HistFac.makeTH1D("AntiEtaZVertex", "z / cm","",bins_ZVertex,"h_AntiEtaZVertex");
        h_TreeZVertex = HistFac.makeTH1D("TreeZVertex", "z / cm","",bins_ZVertex,"h_TreeZVertex");

        h_Pi0EtaFitProb = HistFac.makeTH2D("Pi0 vs. Eta","Pi0 p","Eta p",bins_FitProb,bins_FitProb,"h_Pi0EtaFitProb");
        h_Pi0TreeFitProb = HistFac.makeTH2D("Pi0 vs. Tree","Pi0 p","Tree p",bins_FitProb,bins_FitProb,"h_Pi0TreeFitProb");
        h_EtaTreeFitProb = HistFac.makeTH2D("Eta vs. Tree","Eta p","Tree p",bins_FitProb,bins_FitProb,"h_EtaTreeFitProb");

        BinSettings bins_BachelorE(100,100,200);
        h_Bachelor_E = HistFac.makeTH1D("E_#gamma in #eta' frame","E_{#gamma} / MeV","",
                                               bins_BachelorE,"h_Bachelor_E");

    }

    void Fill(const Fill_t& f) const {
        CommonHist_t::Fill(f);
        const SharedTree_t& s = f.Shared;
        const Tree_t& tree = f.Tree;

        h_IM_4g->Fill(tree.IM_Pi0gg, f.TaggW());

        h_AntiPi0FitProb->Fill(f.AntiPi0FitProb, f.TaggW());
        h_AntiEtaFitProb->Fill(f.AntiEtaFitProb, f.TaggW());
        h_TreeFitProb->Fill(f.TreeFitProb, f.TaggW());

        h_AntiPi0ZVertex->Fill(s.AntiPi0FitZVertex, f.TaggW());
        h_AntiEtaZVertex->Fill(s.AntiEtaFitZVertex, f.TaggW());
        h_TreeZVertex->Fill(tree.TreeFitZVertex, f.TaggW());

        h_Pi0EtaFitProb->Fill(f.AntiPi0FitProb, f.AntiEtaFitProb, f.TaggW());
        h_Pi0TreeFitProb->Fill(f.AntiPi0FitProb, f.TreeFitProb, f.TaggW());
        h_EtaTreeFitProb->Fill(f.AntiEtaFitProb, f.TreeFitProb, f.TaggW());

    }

    std::vector<TH1*> GetHists() const {
        auto hists = CommonHist_t::GetHists();
        hists.insert(hists.end(), {
                         h_IM_4g,
                         h_AntiPi0FitProb, h_AntiEtaFitProb, h_TreeFitProb,
                         h_AntiPi0ZVertex, h_AntiEtaZVertex, h_TreeZVertex,
                         h_Bachelor_E
                     });
        return hists;
    }

    static cuttree::Cuts_t<Fill_t> GetCuts() {
        using cuttree::MultiCut_t;
        auto cuts = cuttree::ConvertCuts<Fill_t, CommonHist_t::Fill_t>(CommonHist_t::GetCuts());

//        cuts.emplace_back(MultiCut_t<Fill_t>{
//                              {"PIDSumE=0", [] (const Fill_t& f) { return f.Common.PIDSumE==0; } },
//                              {"CBSumVetoE=0", [] (const Fill_t& f) { return f.Shared.CBSumVetoE==0; } },
//                              {"NoCut", [] (const Fill_t&) { return true; } },
//                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"AntiPi0FitProb<0.0001||nan", [] (const Fill_t& f)  { return std::isnan(f.Shared.AntiPi0FitProb) || f.Shared.AntiPi0FitProb<0.0001; } },
                              {"AntiPi0FitProb<0.00001||nan", [] (const Fill_t& f) { return std::isnan(f.Shared.AntiPi0FitProb) || f.Shared.AntiPi0FitProb<0.00001; } },
                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"AntiEtaFitProb<0.0001||nan", [] (const Fill_t& f)  { return std::isnan(f.Shared.AntiEtaFitProb) || f.Shared.AntiEtaFitProb<0.0001; } },
                              {"AntiEtaFitProb<0.00001||nan", [] (const Fill_t& f) { return std::isnan(f.Shared.AntiEtaFitProb) || f.Shared.AntiEtaFitProb<0.00001; } },
                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"TreeFitProb>0.2", [] (const Fill_t& f) { return f.Tree.TreeFitProb>0.2; } },
                              {"TreeFitProb>0.1", [] (const Fill_t& f) { return f.Tree.TreeFitProb>0.1; } },
                          });

        auto LowE_cut =  [] (const Fill_t& f)  {
            if(f.Tree.gNonPi0_CaloE().front()<100 || f.Tree.gNonPi0_CaloE().back()<100)
                return false;
            return true;
        };

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"LowE", LowE_cut},
                          });

        auto PID_cut =  [] (const Fill_t& f)  {
            // require both in CB
            if(!f.Tree.gNonPi0_IsCB().front() || !f.Tree.gNonPi0_IsCB().back())
                return false;
            if(f.Tree.gNonPi0_VetoE().front()>1.0 || f.Tree.gNonPi0_VetoE().back()>1.0)
                return false;
            return true;
        };

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"PID", PID_cut},
                          });



        return cuts;
    }
};

struct SigPi0Hist_t : SigHist_t {
    using Tree_t = physics::EtapOmegaG::Sig_t::Pi0_t::BaseTree_t;

    TH2D* h_IM_3g_4g_high;    // Omega IM vs. EtaPrime IM

    struct Fill_t : SigHist_t::Fill_t {
        const Tree_t& Pi0;
        Fill_t(const CommonHist_t::Tree_t& common,
               const SharedTree_t& shared,
               const Tree_t& pi0) :
            SigHist_t::Fill_t(common, shared, pi0),
            Pi0(pi0)
        {}
    };

    SigPi0Hist_t(HistogramFactory HistFac, cuttree::TreeInfo_t treeInfo) : SigHist_t(HistFac, treeInfo) {
        h_IM_3g_4g_high = HistFac.makeTH2D("Best #omega vs. #eta' IM",
                                           "IM(#pi^{0}#gamma#gamma) / MeV",
                                           "IM(#pi^{0}#gamma) / MeV",
                                           bins_IM_Etap, bins_IM_Omega,"h_IM_3g_4g_high_best"
                                           );
    }

    void Fill(const Fill_t& f) const {
        SigHist_t::Fill(f);
        const Tree_t& pi0 = f.Pi0;
        h_IM_3g_4g_high->Fill(pi0.IM_Pi0gg, pi0.IM_Pi0g()[1], f.TaggW());
        h_Bachelor_E->Fill(pi0.Bachelor_E()[0], f.TaggW());
    }

    static cuttree::Cuts_t<Fill_t> GetCuts() {
        using cuttree::MultiCut_t;
        auto cuts = cuttree::ConvertCuts<Fill_t, SigHist_t::Fill_t>(SigHist_t::GetCuts());
        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"IM_Pi0g[1]", [] (const Fill_t& f) {
                                   const auto& window = ParticleTypeDatabase::Omega.GetWindow(80);
                                   return window.Contains(f.Pi0.IM_Pi0g()[1]);
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




    SigOmegaPi0Hist_t(HistogramFactory HistFac, cuttree::TreeInfo_t treeInfo) : SigHist_t(HistFac, treeInfo) {

    }

    std::vector<TH1*> GetHists() const {
        return SigHist_t::GetHists();
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
            CommonHist_t::Fill_t(common, tree),
            Tree(tree)
        {}
    };

    TH1D* h_IM_2g;

    TH2D* h_IM_2g_TaggCh;

    TH1D* h_TaggT;

    RefHist_t(HistogramFactory HistFac, cuttree::TreeInfo_t treeInfo) : CommonHist_t(HistFac, treeInfo) {
        BinSettings bins_im(150,800,1050);

        h_IM_2g = HistFac.makeTH1D("IM 2g","IM / MeV","",bins_im,"h_IM_2g");

        auto ept = ExpConfig::Setup::GetDetector<expconfig::detector::EPT>();
        h_IM_2g_TaggCh = HistFac.makeTH2D("IM 2g vs. TaggCh","IM / MeV","Tagger Channel",
                                          bins_im, BinSettings(ept->GetNChannels()), "h_IM_2g_TaggCh");

        h_TaggT = HistFac.makeTH1D("Tagger Time","t / ns","",BinSettings(400,-50,50),"h_TaggT");
    }

    void Fill(const Fill_t& f) const {
        CommonHist_t::Fill(f);
        const Tree_t& tree = f.Tree;
        h_IM_2g->Fill(tree.IM_2g, f.TaggW());

        h_IM_2g_TaggCh->Fill(tree.IM_2g, f.Common.TaggCh, f.TaggW());

        h_TaggT->Fill(f.Common.TaggT-f.Common.CBAvgTime, f.TaggW()<0 ? -1.0 : 1.0);
    }

    std::vector<TH1*> GetHists() const {
        auto hists = CommonHist_t::GetHists();
        hists.insert(hists.end(), {h_IM_2g});
        return hists;
    }

    static cuttree::Cuts_t<Fill_t> GetCuts() {
        using cuttree::MultiCut_t;
        auto cuts = cuttree::ConvertCuts<Fill_t, CommonHist_t::Fill_t>(CommonHist_t::GetCuts());

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"KinFitProb>0.02", [] (const Fill_t& f) { return f.Shared.KinFitProb>0.02; } },
                              {"KinFitProb>0.05", [] (const Fill_t& f) { return f.Shared.KinFitProb>0.05; } },
                          });
        return cuts;
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

    auto cmd_tree = cmd.add<TCLAP::ValueArg<string>>("t","tree","Tree name",false,"SigOmegaPi0","treename");
    auto cmd_setupname = cmd.add<TCLAP::ValueArg<string>>("s","setup","Override setup name", false, "Setup_2014_07_EPT_Prod", "setup");

    cmd.parse(argc, argv);

    const auto setup_name = cmd_setupname->getValue() ;
    auto setup = ExpConfig::Setup::Get(setup_name);
    if(setup == nullptr) {
        LOG(ERROR) << "Did not find setup instance for name " << setup_name;
        return 1;
    }

    WrapTFileInput input(cmd_input->getValue());

    auto link_branches = [&input] (const string treename, WrapTTree& wraptree, long long expected_entries) {
        TTree* t;
        if(!input.GetObject(treename,t))
            throw runtime_error("Cannot find tree "+treename+" in input file");
        if(expected_entries>=0 && t->GetEntries() != expected_entries)
            throw runtime_error("Tree "+treename+" does not have entries=="+to_string(expected_entries));
        if(wraptree.Matches(t, true, true)) {
            wraptree.LinkBranches(t);
            return true;
        }
        return false;
    };


    CommonHist_t::Tree_t treeCommon;
    if(!link_branches("EtapOmegaG/treeCommon", treeCommon, -1)) {
        LOG(ERROR) << "Cannot link branches of treeCommon";
        return 1;
    }

    auto entries = treeCommon.Tree->GetEntries();



    SigPi0Hist_t::Tree_t treeSigPi0;
    SigOmegaPi0Hist_t::Tree_t treeSigOmegaPi0;
    RefHist_t::Tree_t treeRef;

    // auto-detect which tree "type" to use
    const auto& treename = cmd_tree->getValue();
    if(link_branches("EtapOmegaG/"+treename, treeSigPi0, entries)) {
        LOG(INFO) << "Identified " << treename << " as signal tree (Pi0)";
    }
    else if(link_branches("EtapOmegaG/"+treename, treeSigOmegaPi0, entries)) {
        LOG(INFO) << "Identified " << treename << " as signal tree (OmegaPi0)";
    }
    else if(link_branches("EtapOmegaG/"+treename, treeRef, entries)) {
        LOG(INFO) << "Identified " << treename << " as reference tree";
    }
    else {
        LOG(ERROR) << "Could not identify " << treename;
        return 1;
    }

    SigHist_t::SharedTree_t treeSigShared;
    if(treeSigPi0 || treeSigOmegaPi0) {
        if(!link_branches("EtapOmegaG/SigShared", treeSigShared, entries)) {
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

    long long entry = 0;
    ProgressCounter::Interval = 3;
    ProgressCounter progress(
                [&entry, entries] (std::chrono::duration<double>) {
        LOG(INFO) << "Processed " << 100.0*entry/entries << " %";
    });

    for(entry=0;entry<max_entries;entry++) {
        if(interrupt)
            break;

        treeCommon.Tree->GetEntry(entry);

        // we handle the Ref/Sig cut here to save some reading work
        if(treeSigPi0 || treeSigOmegaPi0) {
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
        else if(treeRef) {
            treeRef.Tree->GetEntry(entry);
            cuttree::Fill<MCRefHist_t>(cuttreeRef, {treeCommon, treeRef});
        }
        ProgressCounter::Tick();
    }

    if(!cmd_batchmode->isSet()) {
        if(!std_ext::system::isInteractive()) {
            LOG(INFO) << "No TTY attached. Not starting ROOT shell.";
        }
        else {

            argc=0; // prevent TRint to parse any cmdline
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
