#include "EtapOmegaG.h"
#include "physics/Plotter.h"

#include "analysis/plot/CutTree.h"

#include "expconfig/ExpConfig.h"
#include "expconfig/detectors/EPT.h"

#include "base/interval.h"
#include "base/Logger.h"

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::plot;
using namespace std;

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
        this->GetHist(5, "D07", Mod_t::MakeDataPoints(kGray));
        this->GetHist(6, "D10", Mod_t::MakeDataPoints(kGray));
        this->GetHist(7, "D12", Mod_t::MakeDataPoints(kGray));

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
                                                      Mod_t::MakeLine(histstyle::color_t::GetLight(mctrue-9), 1, kGray+2)
                                                      );

        hist.Fill(f);

        // handle MC_all and MC_bkg
        if(mctrue>0) {
            this->GetHist(3).Fill(f);
            if(mctrue >= 9)
                this->GetHist(4).Fill(f);
        }

        // handle D07/D10/D12
        if(mctrue == 0) {
            this->GetHist(4+f.Common.BeamTime).Fill(f);
        }
    }
};

// define the structs containing the histograms
// and the cuts. for simple branch variables, that could
// be combined...

struct CommonHist_t {
    static OptionsPtr opts;

    using Tree_t = physics::EtapOmegaG::TreeCommon;
    using ProtonPhoton_t = physics::EtapOmegaG::ProtonPhotonTree_t;
    struct Fill_t {
        const Tree_t& Common;
        const ProtonPhoton_t& ProtonPhoton;
        const utils::MCWeighting::tree_t& MCWeighting;

        Fill_t(const Tree_t& common, const ProtonPhoton_t& protonphoton,
               const utils::MCWeighting::tree_t& mcWeighting) :
            Common(common), ProtonPhoton(protonphoton), MCWeighting(mcWeighting) {}

        double Weight() const {
            if(MCWeighting.Tree)
                return MCWeighting.MCWeight;
            return Common.TaggW;
        }

    };

    const BinSettings bins_FitProb{100, 0, 1};
    const BinSettings bins_LogFitProb{100, -8, -2};
    TH1D* h_CBSumE = nullptr;
    TH1D* h_CBSumVetoE = nullptr;
    TH1D* h_PIDSumE = nullptr;
    TH1D* h_MissingMass = nullptr;
    TH1D* h_DiscardedEk = nullptr;
    TH1D* h_nTouchesHole = nullptr;
    TH1D* h_MCMissedBkg = nullptr;

    TH2D* h_ProtonTOF = nullptr;
    TH2D* h_ProtonTOFFitted = nullptr;
    TH2D* h_ProtonVetoE = nullptr;
    TH2D* h_ProtonShortE = nullptr;

    const bool isLeaf;
    const bool includeProtonHists;


    CommonHist_t(HistogramFactory HistFac, cuttree::TreeInfo_t treeInfo) :
        isLeaf(treeInfo.nDaughters==0),
        includeProtonHists(opts->Get<bool>("IncludeProtonHists", false))
    {
        h_CBSumE = HistFac.makeTH1D("CB Sum E","E / MeV","",BinSettings(100,500,1600),"h_CBSumE");
        h_CBSumVetoE = HistFac.makeTH1D("CB Veto Sum E","E / MeV","",BinSettings(100,0,4),"h_CBSumVetoE");
        h_PIDSumE = HistFac.makeTH1D("PID Sum E","E / MeV","",BinSettings(50,0,10),"h_PIDSumE");
        h_MissingMass = HistFac.makeTH1D("MissingMass","m / MeV","",BinSettings(200,600,1300),"h_MissingMass");
        h_DiscardedEk = HistFac.makeTH1D("DiscardedEk","E / MeV","",BinSettings(100,0,100),"h_DiscardedEk");
        h_nTouchesHole = HistFac.makeTH1D("nTouchesHole","nTouchesHole","",BinSettings(5),"h_nTouchesHole");
        h_MCMissedBkg = HistFac.makeTH1D("MCMissedBkg","","",BinSettings(15),"h_MCMissedBkg");

        if(!isLeaf)
            return;

        if(includeProtonHists) {
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
    }


    void Fill(const Fill_t& f) const {

        h_CBSumE->Fill(f.Common.CBSumE, f.Weight());
        h_CBSumVetoE->Fill(f.ProtonPhoton.CBSumVetoE, f.Weight());
        h_PIDSumE->Fill(f.Common.PIDSumE, f.Weight());
        h_MissingMass->Fill(f.ProtonPhoton.MissingMass, f.Weight());
        h_DiscardedEk->Fill(f.ProtonPhoton.DiscardedEk, f.Weight());
        h_nTouchesHole->Fill(f.ProtonPhoton.nTouchesHole, f.Weight());

        if(!f.Common.MCTrueMissed().empty())
            h_MCMissedBkg->Fill(f.Common.MCTrueMissed().c_str(), f.Weight());

        if(!isLeaf)
            return;

        if(includeProtonHists) {
            h_ProtonTOF->Fill(f.ProtonPhoton.ProtonTime, f.ProtonPhoton.ProtonE, f.Weight());
            h_ProtonTOFFitted->Fill(f.ProtonPhoton.ProtonTime, f.ProtonPhoton.FittedProtonE, f.Weight());
            h_ProtonVetoE->Fill(f.ProtonPhoton.FittedProtonE, f.ProtonPhoton.ProtonVetoE, f.Weight());
            h_ProtonShortE->Fill(f.ProtonPhoton.FittedProtonE, f.ProtonPhoton.ProtonShortE, f.Weight());
        }
    }

    std::vector<TH1*> GetHists() const {
        return {h_CBSumE, h_CBSumVetoE, h_PIDSumE, h_MissingMass,
                    h_DiscardedEk, h_nTouchesHole};
    }

    // Sig and Ref channel (can) share some cuts...
    static cuttree::Cuts_t<Fill_t> GetCuts() {
        using cuttree::MultiCut_t;
        cuttree::Cuts_t<Fill_t> cuts;
        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"DiscardedEk=0", [] (const Fill_t& f) { return f.ProtonPhoton.DiscardedEk == 0; } },
                              {"DiscardedEk<50", [] (const Fill_t& f) { return f.ProtonPhoton.DiscardedEk < 50; } },
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
               const Tree_t& tree,
               const utils::MCWeighting::tree_t& mcWeighting) :
            CommonHist_t::Fill_t(common, tree, mcWeighting),
            Shared(shared),
            Tree(tree)
        {}

    };

    TH1D* h_IM_4g;        // EtaPrime IM

    TH1D* h_KinFitProb;
    TH1D* h_AntiPi0FitProb;
    TH1D* h_AntiEtaFitProb;
    TH1D* h_TreeFitProb;

    TH1D* h_AntiPi0ZVertex;
    TH1D* h_AntiEtaZVertex;
    TH1D* h_TreeZVertex;

    TH2D* h_gNonPi0_CaloE_Theta;
    TH1D* h_gNonPi0_TouchesHoles;
    TH1D* h_gNonPi0_CBSumVetoE;

    TH1D* h_Bachelor_E;

    const BinSettings bins_IM_Etap {100, 910, 1020};
    const BinSettings bins_IM_Omega{100, 700, 900};
    const BinSettings bins_ZVertex{100, -15, 15};

    SigHist_t(HistogramFactory HistFac, cuttree::TreeInfo_t treeInfo) : CommonHist_t(HistFac, treeInfo) {
        h_IM_4g = HistFac.makeTH1D("#eta' IM", "IM(#pi^{0}#gamma#gamma) / MeV","",bins_IM_Etap,"h_IM_4g");

        h_KinFitProb = HistFac.makeTH1D("KinFitProb","p","",bins_FitProb,"h_KinFitProb");
        h_AntiPi0FitProb = HistFac.makeTH1D("AntiPi0FitProb", "log_{10} p","",bins_LogFitProb,"h_AntiPi0FitProb");
        h_AntiEtaFitProb = HistFac.makeTH1D("AntiEtaFitProb", "log_{10} p","",bins_LogFitProb,"h_AntiEtaFitProb");
        h_TreeFitProb = HistFac.makeTH1D("TreeFitProb", "p","",bins_FitProb,"h_TreeFitProb");

        h_AntiPi0ZVertex = HistFac.makeTH1D("AntiPi0ZVertex", "z / cm","",bins_ZVertex,"h_AntiPi0ZVertex");
        h_AntiEtaZVertex = HistFac.makeTH1D("AntiEtaZVertex", "z / cm","",bins_ZVertex,"h_AntiEtaZVertex");
        h_TreeZVertex = HistFac.makeTH1D("TreeZVertex", "z / cm","",bins_ZVertex,"h_TreeZVertex");

        h_gNonPi0_CaloE_Theta = HistFac.makeTH2D("gNonPi0 E_{k} #theta","#theta / #circ","E_{k} / MeV",
                                                 BinSettings(180,0,180), BinSettings(50,0,400), "h_gNonPi0_CaloE_Theta");

        h_gNonPi0_TouchesHoles = HistFac.makeTH1D("gNonPi0_TouchesHole","nTouchesHole","",
                                                  BinSettings(3),"h_gNonPi0_TouchesHole");
        h_gNonPi0_CBSumVetoE = HistFac.makeTH1D("CB Veto Sum E Bachelor Photons","E / MeV","",
                                                BinSettings(100,0,2),"h_gNonPi0_CBSumVetoE");

        BinSettings bins_BachelorE(100,100,200);
        h_Bachelor_E = HistFac.makeTH1D("E_#gamma in #eta' frame","E_{#gamma} / MeV","",
                                        bins_BachelorE,"h_Bachelor_E");


    }

    void Fill(const Fill_t& f) const {
        CommonHist_t::Fill(f);
        const SharedTree_t& s = f.Shared;
        const Tree_t& tree = f.Tree;

        h_IM_4g->Fill(tree.IM_Pi0gg, f.Weight());

        h_KinFitProb->Fill(s.KinFitProb, f.Weight());
        h_AntiPi0FitProb->Fill(std::log10(s.AntiPi0FitProb), f.Weight());
        h_AntiEtaFitProb->Fill(std::log10(s.AntiEtaFitProb), f.Weight());
        h_TreeFitProb->Fill(tree.TreeFitProb, f.Weight());

        h_AntiPi0ZVertex->Fill(s.AntiPi0FitZVertex, f.Weight());
        h_AntiEtaZVertex->Fill(s.AntiEtaFitZVertex, f.Weight());
        h_TreeZVertex->Fill(tree.TreeFitZVertex, f.Weight());

        {
            /// \todo actually the physics class should have converted this... fix this in possible next round
            const auto& theta0 = std_ext::radian_to_degree(tree.gNonPi0_Theta()[0]);
            h_gNonPi0_CaloE_Theta->Fill(theta0, tree.gNonPi0_CaloE()[0], f.Weight());\
            const auto& theta1 = std_ext::radian_to_degree(tree.gNonPi0_Theta()[1]);
            h_gNonPi0_CaloE_Theta->Fill(theta1, tree.gNonPi0_CaloE()[1], f.Weight());
        }

        h_gNonPi0_TouchesHoles->Fill(tree.gNonPi0_TouchesHole()[0]+tree.gNonPi0_TouchesHole()[1], f.Weight());
        h_gNonPi0_CBSumVetoE->Fill(tree.gNonPi0_VetoE()[0]+tree.gNonPi0_VetoE()[1], f.Weight());
    }

    std::vector<TH1*> GetHists() const {
        auto hists = CommonHist_t::GetHists();
        hists.insert(hists.end(), {
                         h_IM_4g, h_KinFitProb,
                         h_AntiPi0FitProb, h_AntiEtaFitProb, h_TreeFitProb,
                         h_AntiPi0ZVertex, h_AntiEtaZVertex, h_TreeZVertex,
                         h_gNonPi0_TouchesHoles, h_gNonPi0_CBSumVetoE,
                         h_Bachelor_E
                     });
        return hists;
    }

    static cuttree::Cuts_t<Fill_t> GetCuts() {
        using cuttree::MultiCut_t;
        auto cuts = cuttree::ConvertCuts<Fill_t, CommonHist_t::Fill_t>(CommonHist_t::GetCuts());


        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"AntiPi0FitProb<0.00001||nan", [] (const Fill_t& f) { return std::isnan(f.Shared.AntiPi0FitProb) || f.Shared.AntiPi0FitProb<0.00001; } },
                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"AntiEtaFitProb<0.0001||nan", [] (const Fill_t& f)  { return std::isnan(f.Shared.AntiEtaFitProb) || f.Shared.AntiEtaFitProb<0.0001; } },
                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"TreeFitProb>0.2", [] (const Fill_t& f) { return f.Tree.TreeFitProb>0.2; } },
                              {"TreeFitProb>0.1", [] (const Fill_t& f) { return f.Tree.TreeFitProb>0.1; } },
                          });

        auto gNonPi0_cut_1 = [] (const Fill_t& f) {
            const auto& theta0 = f.Tree.gNonPi0_Theta()[0];
            const auto& caloE0 = f.Tree.gNonPi0_CaloE()[0];
            const auto& theta1 = f.Tree.gNonPi0_Theta()[1];
            const auto& caloE1 = f.Tree.gNonPi0_CaloE()[1];
            const auto cut = [] (double caloE, double theta) {
                /// \todo fix this in physics class...
                theta = std_ext::radian_to_degree(theta);
                return caloE > 230.0*(1.0-theta/160.0);
            };
            return cut(caloE0, theta0) && cut(caloE1, theta1);
        };

        auto gNonPi0_cut_2 = [] (const Fill_t& f) {
            const auto& theta0 = f.Tree.gNonPi0_Theta()[0];
            const auto& caloE0 = f.Tree.gNonPi0_CaloE()[0];
            const auto& theta1 = f.Tree.gNonPi0_Theta()[1];
            const auto& caloE1 = f.Tree.gNonPi0_CaloE()[1];
            const auto cut = [] (double caloE, double theta) {
                /// \todo fix this in physics class...
                theta = std_ext::radian_to_degree(theta);
                if(theta<22) // decide if TAPS or CB
                    return caloE > 140;
                else
                    return caloE > 60;
            };
            return cut(caloE0, theta0) && cut(caloE1, theta1);
        };

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"gNonPi0_1", gNonPi0_cut_1},
                              {"gNonPi0_2", gNonPi0_cut_2},
                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"CBSumVetoE<0.2", [] (const Fill_t& f) { return f.ProtonPhoton.CBSumVetoE<0.2; }},
                              {"CBSumVetoE_gNonPi0<0.2", [] (const Fill_t& f) {
                                   auto& v = f.Tree.gNonPi0_VetoE();
                                   return (v.front()+v.back())<0.2;
                               }
                              },
                              {"-", [] (const Fill_t&) { return true; }},
                          });
        return cuts;
    }
};

struct SigPi0Hist_t : SigHist_t {
    using Tree_t = physics::EtapOmegaG::Sig_t::Pi0_t::Tree_t;

    TH2D* h_IM_3g_4g_high;    // Omega IM vs. EtaPrime IM

    struct Fill_t : SigHist_t::Fill_t {
        const Tree_t& Pi0;
        Fill_t(const CommonHist_t::Tree_t& common,
               const SharedTree_t& shared,
               const Tree_t& pi0,
               const utils::MCWeighting::tree_t& mcWeighting) :
            SigHist_t::Fill_t(common, shared, pi0, mcWeighting),
            Pi0(pi0)
        {}
    };

    SigPi0Hist_t(HistogramFactory HistFac, cuttree::TreeInfo_t treeInfo) : SigHist_t(HistFac, treeInfo) {
        h_IM_3g_4g_high = HistFac.makeTH2D("#omega vs. #eta' IM",
                                           "IM(#pi^{0}#gamma#gamma) / MeV",
                                           "IM(#pi^{0}#gamma) / MeV",
                                           bins_IM_Etap, bins_IM_Omega,"h_IM_3g_4g_high"
                                           );
    }

    void Fill(const Fill_t& f) const {
        SigHist_t::Fill(f);
        const Tree_t& pi0 = f.Pi0;
        h_IM_3g_4g_high->Fill(pi0.IM_Pi0gg, pi0.IM_Pi0g()[1], f.Weight());
        h_Bachelor_E->Fill(pi0.Bachelor_E()[0], f.Weight());
    }

    static cuttree::Cuts_t<Fill_t> GetCuts() {
        using cuttree::MultiCut_t;
        auto cuts = cuttree::ConvertCuts<Fill_t, SigHist_t::Fill_t>(SigHist_t::GetCuts());
        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"IM_Pi0g[1]", [] (const Fill_t& f) {
                                   const auto& window = ParticleTypeDatabase::Omega.GetWindow(40);
                                   return window.Contains(f.Pi0.IM_Pi0g()[1]);
                               }},
                          });
        return cuts;
    }

};

struct SigOmegaPi0Hist_t : SigHist_t {
    using Tree_t = physics::EtapOmegaG::Sig_t::OmegaPi0_t::Tree_t;
    struct Fill_t : SigHist_t::Fill_t {
        const Tree_t& OmegaPi0;
        Fill_t(const CommonHist_t::Tree_t& common,
               const SharedTree_t& shared,
               const Tree_t& omegapi0,
               const utils::MCWeighting::tree_t& mcWeighting) :
            SigHist_t::Fill_t(common, shared, omegapi0, mcWeighting),
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
        h_Bachelor_E->Fill(omegapi0.Bachelor_E, f.Weight());
    }

    static cuttree::Cuts_t<Fill_t> GetCuts() {
        return cuttree::ConvertCuts<Fill_t, SigHist_t::Fill_t>(SigHist_t::GetCuts());
    }

};

struct RefHist_t : CommonHist_t {
    using Tree_t = physics::EtapOmegaG::Ref_t::Tree_t;

    struct Fill_t : CommonHist_t::Fill_t {
        const Tree_t& Tree;
        Fill_t(const CommonHist_t::Tree_t& common, const Tree_t& tree,
               const utils::MCWeighting::tree_t& mcWeighting) :
            CommonHist_t::Fill_t(common, tree, mcWeighting),
            Tree(tree)
        {}

    };


    TH1D* h_KinFitProb;
    TH1D* h_IM_2g;
    TH2D* h_IM_2g_TaggCh;

    TH1D* h_TaggT;

    RefHist_t(HistogramFactory HistFac, cuttree::TreeInfo_t treeInfo) : CommonHist_t(HistFac, treeInfo) {
        BinSettings bins_im(150,875,1050);

        h_KinFitProb = HistFac.makeTH1D("KinFitProb","p","",bins_FitProb,"h_KinFitProb");
        h_IM_2g = HistFac.makeTH1D("IM 2g","IM / MeV","",bins_im,"h_IM_2g");

        auto ept = ExpConfig::Setup::GetDetector<expconfig::detector::EPT>();
        h_IM_2g_TaggCh = HistFac.makeTH2D("IM 2g vs. TaggCh","IM / MeV","Tagger Channel",
                                          bins_im, BinSettings(ept->GetNChannels()), "h_IM_2g_TaggCh");

        h_TaggT = HistFac.makeTH1D("Tagger Time","t / ns","",BinSettings(400,-50,50),"h_TaggT");
    }

    void Fill(const Fill_t& f) const {
        CommonHist_t::Fill(f);
        const Tree_t& tree = f.Tree;

        h_KinFitProb->Fill(tree.KinFitProb, f.Weight());

        h_IM_2g->Fill(tree.IM_2g, f.Weight());
        h_IM_2g_TaggCh->Fill(tree.IM_2g, f.Common.TaggCh, f.Weight());

        h_TaggT->Fill(f.Common.TaggT-f.Common.CBAvgTime, f.Weight()<0 ? -1.0 : 1.0);
    }

    std::vector<TH1*> GetHists() const {
        auto hists = CommonHist_t::GetHists();
        hists.insert(hists.end(), {h_KinFitProb, h_IM_2g});
        return hists;
    }

    static cuttree::Cuts_t<Fill_t> GetCuts() {
        using cuttree::MultiCut_t;
        auto cuts = cuttree::ConvertCuts<Fill_t, CommonHist_t::Fill_t>(CommonHist_t::GetCuts());

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"KinFitProb>0.02", [] (const Fill_t& f) { return f.Tree.KinFitProb>0.02; } },
                              {"KinFitProb>0.05", [] (const Fill_t& f) { return f.Tree.KinFitProb>0.05; } },
                          });
        return cuts;
    }
};

OptionsPtr CommonHist_t::opts;

struct EtapOmegaG_plot : Plotter {

    CommonHist_t::Tree_t treeCommon;
    utils::MCWeighting::tree_t treeMCWeighting;

    EtapOmegaG_plot(const string& tag, const string& name, const WrapTFileInput& input, OptionsPtr opts) :
        Plotter(name, input, opts)
    {
        /// \todo using a static field is actually quite ugly...
        CommonHist_t::opts = opts;

        init_tree(input, treeCommon, "EtapOmegaG/"+tag+"/Common");

        if(input.GetObject("EtapOmegaG/"+tag+"/"+utils::MCWeighting::treeName, treeMCWeighting.Tree)) {
            LOG(INFO) << "Found " << tag << " MCWeighting tree";
            treeMCWeighting.LinkBranches();
            check_entries(treeMCWeighting);
        }
    }

    static void init_tree(const WrapTFileInput& input, WrapTTree& tree, const string& name) {
        if(!input.GetObject(name, tree.Tree))
            throw Exception("Cannot find tree "+name);
        tree.LinkBranches();
    }

    void check_entries(const WrapTTree& tree) {
        if(tree.Tree->GetEntries() == treeCommon.Tree->GetEntries())
            return;
        throw Exception(std_ext::formatter() << "Tree " << tree.Tree->GetName()
                        << " does not have expected entries=" << treeCommon.Tree->GetEntries()
                        << " but " << tree.Tree->GetEntries());
    }

    virtual long long GetNumEntries() const override
    {
        return treeCommon.Tree->GetEntries();
    }

    virtual void ProcessEntry(const long long entry) override
    {
        treeCommon.Tree->GetEntry(entry);
        if(treeMCWeighting.Tree)
            treeMCWeighting.Tree->GetEntry(entry);
    }

};

struct EtapOmegaG_plot_Ref : EtapOmegaG_plot {

    RefHist_t::Tree_t          treeRef;

    using MCRefHist_t = MCTrue_Splitter<RefHist_t>;
    cuttree::Tree_t<MCRefHist_t> cuttreeRef;

    EtapOmegaG_plot_Ref(const string& name, const WrapTFileInput& input, OptionsPtr opts) :
        EtapOmegaG_plot("Ref", name, input, opts)
    {
        init_tree(input, treeRef, "EtapOmegaG/Ref/Ref");
        check_entries(treeRef);

        cuttreeRef = cuttree::Make<MCRefHist_t>(HistFac);
    }

    virtual void ProcessEntry(const long long entry) override
    {
        EtapOmegaG_plot::ProcessEntry(entry);
        treeRef.Tree->GetEntry(entry);
        cuttree::Fill<MCRefHist_t>(cuttreeRef, {treeCommon, treeRef, treeMCWeighting});
    }
};

struct EtapOmegaG_plot_Sig : EtapOmegaG_plot {

    SigHist_t::SharedTree_t   treeSigShared;
    SigPi0Hist_t::Tree_t      treeSigPi0;
    SigOmegaPi0Hist_t::Tree_t treeSigOmegaPi0;

    using MCSigPi0Hist_t = MCTrue_Splitter<SigPi0Hist_t>;
    using MCSigOmegaPi0Hist_t = MCTrue_Splitter<SigOmegaPi0Hist_t>;
    cuttree::Tree_t<MCSigPi0Hist_t>      cuttreeSigPi0;
    cuttree::Tree_t<MCSigOmegaPi0Hist_t> cuttreeSigOmegaPi0;

    EtapOmegaG_plot_Sig(const string& name, const WrapTFileInput& input, OptionsPtr opts) :
        EtapOmegaG_plot("Sig", name, input, opts)
    {
        init_tree(input, treeSigShared,   "EtapOmegaG/Sig/Shared");
        init_tree(input, treeSigPi0,      "EtapOmegaG/Sig/Pi0");
        init_tree(input, treeSigOmegaPi0, "EtapOmegaG/Sig/OmegaPi0");

        check_entries(treeSigShared);
        check_entries(treeSigPi0);
        check_entries(treeSigOmegaPi0);

        cuttreeSigPi0 = cuttree::Make<MCSigPi0Hist_t>(HistogramFactory("SigPi0",HistFac,"SigPi0"));
        cuttreeSigOmegaPi0 = cuttree::Make<MCSigOmegaPi0Hist_t>(HistogramFactory("SigOmegaPi0",HistFac,"SigOmegaPi0"));
    }

    virtual void ProcessEntry(const long long entry) override
    {
        EtapOmegaG_plot::ProcessEntry(entry);
        treeSigShared.Tree->GetEntry(entry);
        treeSigPi0.Tree->GetEntry(entry);
        treeSigOmegaPi0.Tree->GetEntry(entry);
        cuttree::Fill<MCSigPi0Hist_t>(cuttreeSigPi0, {treeCommon, treeSigShared, treeSigPi0, treeMCWeighting});
        cuttree::Fill<MCSigOmegaPi0Hist_t>(cuttreeSigOmegaPi0, {treeCommon, treeSigShared, treeSigOmegaPi0, treeMCWeighting});
    }
};

AUTO_REGISTER_PLOTTER(EtapOmegaG_plot_Ref)
AUTO_REGISTER_PLOTTER(EtapOmegaG_plot_Sig)
