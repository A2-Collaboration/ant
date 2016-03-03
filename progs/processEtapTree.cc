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
    SetupLogger();

    TCLAP::CmdLine cmd("process analysis tree for physics class Etap3Pi0", ' ', "0.1");


    auto cmd_infile = cmd.add<TCLAP::ValueArg<string>>("i", "i", "File containing tree", false, "", "filename");

    cmd.parse(argc, argv);


    int fake_argc=1;
    char* fake_argv[2];
    fake_argv[0] = argv[0];
    auto app = new TRint("Ant-calib",&fake_argc,fake_argv);


    HistogramFactory HistFac("hists");

    BinSettings singlePion(250,0,700);
    BinSettings doublePion(360,0,1000);

    TH2D* all      = HistFac.makeTH2D("All events",                   "m(2 #pi^{0}) [MeV]","m(#pi^{0}) [MeV]",doublePion,singlePion,"all");
    TH2D* signal   = HistFac.makeTH2D("Signal selection",             "m(2 #pi^{0}) [MeV]","m(#pi^{0}) [MeV]",doublePion,singlePion,"signal");
    TH2D* ref      = HistFac.makeTH2D("Reference selection",          "m(2 #pi^{0}) [MeV]","m(#pi^{0}) [MeV]",doublePion,singlePion,"ref");
    TH2D* mcSignal = HistFac.makeTH2D("Signal mc true",               "m(2 #pi^{0}) [MeV]","m(#pi^{0}) [MeV]",doublePion,singlePion,"mcSignal");
    TH2D* mcRef    = HistFac.makeTH2D("Reference mc true",            "m(2 #pi^{0}) [MeV]","m(#pi^{0}) [MeV]",doublePion,singlePion,"mcRef");
    TH2D* mc3pi0   = HistFac.makeTH2D("3 #pi^{0} production mc true", "m(2 #pi^{0}) [MeV]","m(#pi^{0}) [MeV]",doublePion,singlePion,"mc3pi0");

    WrapTFileInput file(cmd_infile->getValue());
    TTree* tree;
    file.GetObject("Etap3pi0/tree",tree);

    double taggW;
    int type;
    int truetype;
    vector<TLorentzVector> pions(3);
    vector<TLorentzVector> gammas(6);
    TLorentzVector sixG;
    TLorentzVector proton;

    tree->SetBranchAddress("taggWeight", &taggW);
    tree->SetBranchAddress("type",&type);
    tree->SetBranchAddress("truetype",&truetype);
    tree->SetBranchAddress("kf_inter_Sig", &pions);
    tree->SetBranchAddress("kf_gammas_Sig", &gammas);
    tree->SetBranchAddress("kf_6g", &sixG);
    tree->SetBranchAddress("kf_p", &proton);

    for (std::int64_t ind_entry = 0 ; ind_entry < tree->GetEntries() ; ++ind_entry)
    {
        tree->GetEntry(ind_entry);
        for ( unsigned i = 0 ; i < 3 ; ++i)
        {
//            const auto pM  = gammas[i]->M();
//            const auto ppM = (*gammas.at((i + 2) % 3) + *gammas.at((i + 1) % 3)).M();
            /*
            all->Fill(ppM, pM, taggW);
            if ( type == 0 )
                signal->Fill(ppM, pM, taggW);
            if ( type == 1 )
                ref->Fill(ppM, pM, taggW);
            if ( truetype == 0 )
                mcSignal->Fill(ppM, pM, taggW);
            if ( truetype == 1 )
                mcRef->Fill(ppM, pM, taggW);
            if ( truetype == 2 )
                mc3pi0->Fill(ppM, pM, taggW);
                */
        }
    }

    gStyle->SetOptStat(kFALSE);
    canvas("identified") << drawoption("col") << all << signal << ref << endc;
    canvas("mctrue") << drawoption("col") << mcSignal << mcRef << mc3pi0 << endc;

    app->Run(kTRUE);

}

