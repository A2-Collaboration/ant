#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include "base/WrapTFile.h"

#include "analysis/plot/HistogramFactory.h"
#include "TH2D.h"
#include "TTree.h"
#include "TRint.h"
#include "TStyle.h"
#include "TLorentzVector.h"


#include "tclap/CmdLine.h"
#include "base/std_ext/string.h"
#include "base/Logger.h"
#include "tclap/ValuesConstraintExtra.h"
#include "base/std_ext/iterators.h"

#include "analysis/plot/RootDraw.h"

#include "root-addons/analysis_codes/hstack.h"

#include "analysis/utils/Combinatorics.h"


using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::utils;

struct tree_reader_t
{
    tree_reader_t(TTree* theTree):
        tree(theTree)
    {}
    TTree* tree;

    double taggW;
    int type;
    int truetype;
    vector<TLorentzVector>* pions = new vector<TLorentzVector>(3);
    vector<TLorentzVector>* gammas = new vector<TLorentzVector>(6);
    TLorentzVector* sixG = new TLorentzVector();
    TLorentzVector* proton = new TLorentzVector();

    void SetBranches()
    {
        if (tree)
        {
            tree->SetBranchAddress("taggWeight", &taggW);
            tree->SetBranchAddress("type",&type);
            tree->SetBranchAddress("truetype",&truetype);
            tree->SetBranchAddress("kf_inter_Sig", &pions);
            tree->SetBranchAddress("kf_gammas_Sig", &gammas);
            tree->SetBranchAddress("kf_6g", &sixG);
            tree->SetBranchAddress("kf_p", &proton);
        }
    }
};

int main( int argc, char** argv )
{
    // ===================== INIT  ==============================================================================
    SetupLogger();

    TCLAP::CmdLine cmd("process analysis tree for physics class Etap3Pi0", ' ', "0.1");


    auto cmd_infile  = cmd.add<TCLAP::ValueArg<string>>("i", "i", "File containing tree", false, "", "filename");
    auto cmd_infile2 = cmd.add<TCLAP::ValueArg<string>>("c", "c", "File containing tree for comparison", false, "", "filename");

    cmd.parse(argc, argv);

    int fake_argc=1;
    char* fake_argv[2];
    fake_argv[0] = argv[0];
    auto app = new TRint("Ant-calib",&fake_argc,fake_argv);



    // ===================== HIST  ==============================================================================
    HistogramFactory HistFac("hists");

//    BinSettings singlePion(360,0,700);
    BinSettings doublePion(360,0,1000);
    BinSettings doublePionSqr(360,0,800000);
    BinSettings PionP(360,900,1900);
    BinSettings doublePionP(360,1000,2000);

    std::map<std::string,std::map<std::string,TH2D*>> hists2d;
    std::map<std::string,std::map<std::string,TH1D*>> hists1d;

    string cat("pionsID");
    hists2d[cat]["all"]      = HistFac.makeTH2D("All events",                   "m^{2}(4 #gamma) [MeV^{2}]","m^{2}(4 #gamma) [MeV^{2}]",doublePionSqr,doublePionSqr,"all");
    hists2d[cat]["signal"]   = HistFac.makeTH2D("Signal selection",             "m^{2}(4 #gamma) [MeV^{2}]","m^{2}(4 #gamma) [MeV^{2}]",doublePionSqr,doublePionSqr,"signal");
    hists2d[cat]["ref"]      = HistFac.makeTH2D("Reference selection",          "m^{2}(4 #gamma) [MeV^{2}]","m^{2}(4 #gamma) [MeV^{2}]",doublePionSqr,doublePionSqr,"ref");
    cat = "pionsMC";
    hists2d[cat]["mcSignal"] = HistFac.makeTH2D("Signal mc true",               "m^{2}(4 #gamma) [MeV]","m^{2}(4 #gamma) [MeV]",doublePionSqr,doublePionSqr,"mcSignal");
    hists2d[cat]["mcRef"]    = HistFac.makeTH2D("Reference mc true",            "m^{2}(4 #gamma) [MeV]","m^{2}(4 #gamma) [MeV]",doublePionSqr,doublePionSqr,"mcRef");
    hists2d[cat]["mc3pi0"]   = HistFac.makeTH2D("3 #pi^{0} production mc true", "m^{2}(4 #gamma) [MeV]","m^{2}(4 #gamma) [MeV]",doublePionSqr,doublePionSqr,"mc3pi0");
    cat = "NucleonRes";
    hists2d[cat]["1pi"]      = HistFac.makeTH2D("Invariant Mass",               "m(p #pi^{0}) [MeV]","m(2 #pi^{0}) [MeV]",PionP,doublePion,"1pi");
    hists1d[cat]["2pi"]      = HistFac.makeTH1D("Invariant Mass",               "m(p 2 #pi^{0}) [MeV]","#",doublePionP,"2pi");

    const vector<pair<size_t,size_t>> combinations = { { 0 , 1 } , { 0 , 2 } , { 1 , 2 } };




    // ===================== PLOT  ==============================================================================
    if (!cmd_infile2->isSet())
    {
        WrapTFileInput file(cmd_infile->getValue());
        TTree* tree;
        file.GetObject("Etap3pi0/tree",tree);
        tree_reader_t vars(tree);
        vars.SetBranches();

        for (std::int64_t ind_entry = 0 ; ind_entry < tree->GetEntries() ; ++ind_entry)
        {
            tree->GetEntry(ind_entry);
            if ( vars.sixG->M() < 850 )
                continue;
            for ( size_t i = 0 ; i < 3 ; ++i)
                for ( size_t j = 0 ; j < 3 ; ++j)
                {
                    if ( i == j )
                        continue;
                    const auto ppM2  =(vars.pions->at(combinations.at(i).first) + vars.pions->at(combinations.at(i).second)).M2();
                    const auto ppM1  =(vars.pions->at(combinations.at(j).first) + vars.pions->at(combinations.at(j).second)).M2();
                    hists2d.at("pionsID").at("all")->Fill(ppM2,ppM1,vars.taggW);
                    if ( vars.type == 0 )
                        hists2d.at("pionsID").at("signal")->Fill(ppM2, ppM1, vars.taggW);
                    if ( vars.type == 1 )
                        hists2d.at("pionsID").at("ref")->Fill(ppM2, ppM1, vars.taggW);
                    if ( vars.truetype == 0 )
                        hists2d.at("pionsMC").at("mcSignal")->Fill(ppM2, ppM1, vars.taggW);
                    if ( vars.truetype == 1 )
                        hists2d.at("pionsMC").at("mcRef")->Fill(ppM2, ppM1, vars.taggW);
                    if ( vars.truetype == 2 )
                        hists2d.at("pionsMC").at("mc3pi0")->Fill(ppM2, ppM1, vars.taggW);
                }
            for(auto comb=utils::makeCombination(*vars.pions, 3); !comb.done(); ++comb)
            {
                const auto N    = comb.at(0) + *(vars.proton);
                const auto pipi = comb.at(1) + comb.at(2);
                hists2d.at("NucleonRes").at("1pi")->Fill(N.M(),pipi.M(),vars.taggW);
            }
            for(auto comb=utils::makeCombination(*vars.pions, 2); !comb.done(); ++comb)
            {
                const auto N  = comb.at(0) + comb.at(1) + *(vars.proton);
                hists1d.at("NucleonRes").at("2pi")->Fill(N.M(),vars.taggW);
            }

        }
        // ===================== DRAW  ==============================================================================
        gStyle->SetOptStat(kFALSE);

        vector<string> cats2d = {"pionsID","pionsMC"};
        for (auto& category: cats2d)
        {
            canvas c(category);
            c << drawoption("colz");
            for (auto h: hists2d.at(category))
                c << h.second;
            c << endc;
        }
        canvas c("NucleonRes");
        c << hists1d.at("NucleonRes").at("2pi") << drawoption("colz") << hists2d.at("NucleonRes").at("1pi") << endc;
    }
    else
    {
        WrapTFileInput file(cmd_infile->getValue());
        TTree* tree;
        file.GetObject("Etap3pi0/tree",tree);
        WrapTFileInput fileC(cmd_infile2->getValue());
        TTree* tree_comp = nullptr;
        fileC.GetObject("Etap3pi0/tree",tree_comp);

        tree_reader_t vars(tree);
        vars.SetBranches();
        tree_reader_t vars_comp(tree_comp);
        vars_comp.SetBranches();

        BinSettings bins6g(400,400,1050);
        cat = "mc";
        hists1d[cat]["SIGall"]      = HistFac.makeTH1D("all events after signal-cuts", "m(6 #gamma) [MeV]","#",bins6g);
        hists1d[cat]["SIGNoSIG"]    = HistFac.makeTH1D("without signal", "m(6 #gamma) [MeV]","#",bins6g);
        hists1d[cat]["SIG3pi0"]        = HistFac.makeTH1D("3#pi^{0} production", "m(6 #gamma) [MeV]","#",bins6g);
        hists1d[cat]["SIGREF"]        = HistFac.makeTH1D("reference", "m(6 #gamma) [MeV]","#",bins6g);
        hists1d[cat]["3pi0"]        = HistFac.makeTH1D("mc true 3#pi^{0} production", "m(6 #gamma) [MeV]","#",bins6g);
        hists1d[cat]["REF"]         = HistFac.makeTH1D("mc true reference", "m(6 #gamma) [MeV]","#",bins6g);
        hists1d[cat]["SIG"]         = HistFac.makeTH1D("mc true signal", "m(6 #gamma) [MeV]","#",bins6g);
        hists1d[cat]["REFall"]      = HistFac.makeTH1D("all events after reference-cuts", "m(6 #gamma) [MeV]","#",bins6g);
        hists1d[cat]["REFNoREF"]    = HistFac.makeTH1D("without reference", "m(6 #gamma) [MeV]","#",bins6g);
        hists1d[cat]["REF3pi0"]    = HistFac.makeTH1D("3#pi^{0} production", "m(6 #gamma) [MeV]","#",bins6g);
        cat = "data";
        hists1d[cat]["SIGsel"]      = HistFac.makeTH1D("data: identified signal", "m(6 #gamma) [MeV]","#",bins6g);
        hists1d[cat]["REFsel"]      = HistFac.makeTH1D("data: identified reference", "m(6 #gamma) [MeV]","#",bins6g);
        for (std::int64_t ind_entry = 0 ; ind_entry < tree->GetEntries() ; ++ind_entry)
        {
            tree->GetEntry(ind_entry);
            const auto tw = vars.taggW;
            if (vars.type==0)
            {
                hists1d.at("mc").at("SIGall")->Fill(vars.sixG->M(),tw);
                if(vars.truetype!=0)
                    hists1d.at("mc").at("SIGNoSIG")->Fill(vars.sixG->M(),tw);
                if(vars.truetype==1)
                    hists1d.at("mc").at("SIGREF")->Fill(vars.sixG->M(),tw);
                if(vars.truetype==2)
                    hists1d.at("mc").at("SIG3pi0")->Fill(vars.sixG->M(),tw);
            }
            if (vars.type==1)
            {
                hists1d.at("mc").at("REFall")->Fill(vars.sixG->M(),tw);
                if(vars.truetype!=1)
                    hists1d.at("mc").at("REFNoREF")->Fill(vars.sixG->M(),tw);
                if(vars.truetype==2)
                    hists1d.at("mc").at("REF3pi0")->Fill(vars.sixG->M(),tw);
            }
            if (vars.truetype == 0)
                hists1d.at("mc").at("SIG")->Fill(vars.sixG->M(),tw);
            if (vars.truetype == 1)
                hists1d.at("mc").at("REF")->Fill(vars.sixG->M(),tw);
            if (vars.truetype == 2)
                hists1d.at("mc").at("3pi0")->Fill(vars.sixG->M(),tw);
        }
        for (std::int64_t ind_entry = 0 ; ind_entry < tree_comp->GetEntries() ; ++ind_entry)
        {
            tree_comp->GetEntry(ind_entry);
            const auto tw = vars_comp.taggW;
            if (vars_comp.type == 0)
                hists1d.at("data").at("SIGsel")->Fill(vars_comp.sixG->M(),tw);
            if (vars_comp.type == 1)
                hists1d.at("data").at("REFsel")->Fill(vars_comp.sixG->M(),tw);
        }
        canvas s("signal");
        hstack ss("stackSignal","");
        ss << hists1d.at("mc").at("SIGall")
           << hists1d.at("mc").at("SIGNoSIG")
           << hists1d.at("mc").at("SIG3pi0")
           << hists1d.at("mc").at("SIGREF")
           << hists1d.at("mc").at("SIG")
           << hists1d.at("data").at("SIGsel");
        s << drawoption("nostack") << padoption::Legend
          << &ss << endc;

        canvas r("reference");
        hstack sr("stackReference","");
        sr << hists1d.at("mc").at("REFall")
           << hists1d.at("mc").at("REFNoREF")
           << hists1d.at("mc").at("REF3pi0")
           << hists1d.at("data").at("REFsel");
        r << drawoption("nostack") << padoption::Legend
          << &sr << endc;
    }



    app->Run(kTRUE);
}

