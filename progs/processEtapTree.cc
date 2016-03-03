#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include "base/WrapTFile.h"

#include "analysis/plot/HistogramFactories.h"
#include "TH2D.h"
#include "TTree.h"
#include "TRint.h"
#include "TStyle.h"
#include "TLorentzVector.h"


#include "base/CmdLine.h"
#include "base/std_ext/string.h"
#include "base/Logger.h"
#include "base/detail/tclap/ValuesConstraintExtra.h"
#include "base/std_ext/iterators.h"

#include "analysis/plot/root_draw.h"

#include "root-addons/analysis_codes/hstack.h"

#include "analysis/utils/combinatorics.h"


using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::utils;

int main( int argc, char** argv )
{
    // ===================== INIT  ==============================================================================
    SetupLogger();

    TCLAP::CmdLine cmd("process analysis tree for physics class Etap3Pi0", ' ', "0.1");


    auto cmd_infile = cmd.add<TCLAP::ValueArg<string>>("i", "i", "File containing tree", false, "", "filename");

    cmd.parse(argc, argv);


    int fake_argc=1;
    char* fake_argv[2];
    fake_argv[0] = argv[0];
    auto app = new TRint("Ant-calib",&fake_argc,fake_argv);



    // ===================== HIST  ==============================================================================
    HistogramFactory HistFac("hists");

    BinSettings singlePion(360,0,700);
    BinSettings doublePion(360,0,1000);
    BinSettings PionP(360,900,1900);
    BinSettings doublePionP(360,1000,2000);

    std::map<std::string,std::map<std::string,TH2D*>> hists2d;

    string cat("pionsID");
    hists2d[cat]["all"]      = HistFac.makeTH2D("All events",                   "m(2 #pi^{0}) [MeV]","m(#pi^{0}) [MeV]",doublePion,singlePion,"all");
    hists2d[cat]["signal"]   = HistFac.makeTH2D("Signal selection",             "m(2 #pi^{0}) [MeV]","m(#pi^{0}) [MeV]",doublePion,singlePion,"signal");
    hists2d[cat]["ref"]      = HistFac.makeTH2D("Reference selection",          "m(2 #pi^{0}) [MeV]","m(#pi^{0}) [MeV]",doublePion,singlePion,"ref");
    cat = "pionsMC";
    hists2d[cat]["mcSignal"] = HistFac.makeTH2D("Signal mc true",               "m(2 #pi^{0}) [MeV]","m(#pi^{0}) [MeV]",doublePion,singlePion,"mcSignal");
    hists2d[cat]["mcRef"]    = HistFac.makeTH2D("Reference mc true",            "m(2 #pi^{0}) [MeV]","m(#pi^{0}) [MeV]",doublePion,singlePion,"mcRef");
    hists2d[cat]["mc3pi0"]   = HistFac.makeTH2D("3 #pi^{0} production mc true", "m(2 #pi^{0}) [MeV]","m(#pi^{0}) [MeV]",doublePion,singlePion,"mc3pi0");
    cat = "NucleonRes";
    hists2d[cat]["1pi"]      = HistFac.makeTH2D("Invariant Mass",               "m(p #pi^{0}) [MeV]","m(2 #pi^{0}) [MeV]",PionP,doublePion,"1pi");
    hists2d[cat]["2pi"]      = HistFac.makeTH2D("Invariant Mass",               "m(p 2 #pi^{0}) [MeV]","m(#pi^{0}) [MeV]",doublePionP,singlePion,"2pi");


    // ===================== TREE  ==============================================================================
    WrapTFileInput file(cmd_infile->getValue());
    TTree* tree;
    file.GetObject("Etap3pi0/tree",tree);

    double taggW;
    int type;
    int truetype;
    vector<TLorentzVector>* pions = new vector<TLorentzVector>(3);
    vector<TLorentzVector>* gammas = new vector<TLorentzVector>(6);
    TLorentzVector* sixG = new TLorentzVector();
    TLorentzVector* proton = new TLorentzVector();

    tree->SetBranchAddress("taggWeight", &taggW);
    tree->SetBranchAddress("type",&type);
    tree->SetBranchAddress("truetype",&truetype);
    tree->SetBranchAddress("kf_inter_Sig", &pions);
    tree->SetBranchAddress("kf_gammas_Sig", &gammas);
    tree->SetBranchAddress("kf_6g", &sixG);
    tree->SetBranchAddress("kf_p", &proton);

    // ===================== PLOT  ==============================================================================
    for (std::int64_t ind_entry = 0 ; ind_entry < tree->GetEntries() ; ++ind_entry)
    {
        tree->GetEntry(ind_entry);
        for ( unsigned i = 0 ; i < 3 ; ++i)
        {
            const auto pM  = pions->at(i).M();
            const auto ppM = (pions->at((i + 2) % 3) + pions->at((i + 1) % 3)).M();
            hists2d.at("pionsID").at("all")->Fill(ppM,pM,taggW);
            if ( type == 0 )
                hists2d.at("pionsID").at("signal")->Fill(ppM, pM, taggW);
            if ( type == 1 )
                hists2d.at("pionsID").at("ref")->Fill(ppM, pM, taggW);
            if ( truetype == 0 )
                hists2d.at("pionsMC").at("mcSignal")->Fill(ppM, pM, taggW);
            if ( truetype == 1 )
                hists2d.at("pionsMC").at("mcRef")->Fill(ppM, pM, taggW);
            if ( truetype == 2 )
                hists2d.at("pionsMC").at("mc3pi0")->Fill(ppM, pM, taggW);
        }
        for(auto comb=utils::makeCombination(*pions, 3); !comb.Done(); ++comb)
        {
            const auto N    = comb.at(0) + *(proton);
            const auto pipi = comb.at(1) + comb.at(2);
            hists2d.at("NucleonRes").at("1pi")->Fill(N.M(),pipi.M(),taggW);
        }
        for(auto comb=utils::makeCombination(*pions, 2); !comb.Done(); ++comb)
        {
            const auto N  = comb.at(0) + comb.at(1) + *(proton);
            const auto pi = *comb.begin_not();
            hists2d.at("NucleonRes").at("2pi")->Fill(N.M(),pi.M(),taggW);
        }
    }



    // ===================== DRAW  ==============================================================================
    gStyle->SetOptStat(kFALSE);

    vector<string> cats2d = {"pionsID","pionsMC","NucleonRes"};
    for (auto& category: cats2d)
    {
        canvas c(category);
        c << drawoption("colz");
        for (auto h: hists2d.at(category))
            c << h.second;
        c << endc;
    }
    app->Run(kTRUE);

}

