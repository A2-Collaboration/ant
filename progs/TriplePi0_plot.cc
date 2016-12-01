#include <string>
#include <map>
#include <algorithm>

#include "analysis/physics/production/triplePi0.h"
#include "analysis/plot/root_draw.h"

#include "base/CmdLine.h"
#include "base/Logger.h"
#include "base/PlotExt.h"
#include "base/std_ext/math.h"
#include "base/std_ext/string.h"
#include "base/std_ext/system.h"
#include "base/WrapTFile.h"
#include "base/TH_ext.h"

#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMultiGraph.h"
#include "TRint.h"
#include "TTree.h"
#include "TGraphErrors.h"



using namespace ant;
using namespace std;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace ant::std_ext;
using namespace ant::TH_ext;



static volatile bool interrupt = false;
static bool histOut = false;


auto failExit = [] (const string& message)
{
    LOG(ERROR) << message;
    exit(EXIT_FAILURE);
};




int main( int argc, char** argv )
{
    SetupLogger();

    signal(SIGINT, [] (int) {
        LOG(INFO) << ">>> Interrupted";
        interrupt = true;
    });

    TCLAP::CmdLine cmd("TriplePi0_plot", ' ', "0.1");

    auto cmd_input      = cmd.add<TCLAP::ValueArg<string>>("i", "input",     "Input file from triplePi0 class",true,"","filename");
    auto cmd_output     = cmd.add<TCLAP::ValueArg<string>>("o", "save-hist", "Save results to a histogram", false, "","filename");

    //switches
    auto cmd_batchmode  = cmd.add<TCLAP::SwitchArg>("b", "batch", "Run in batch mode (no ROOT shell afterwards)");

    cmd.parse(argc, argv);


    WrapTFileInput inFile(cmd_input->getValue());
    auto treeName = triplePi0::treeAccessName();
    triplePi0::PionProdTree tree;
    if (!inFile.GetObject(treeName,tree.Tree))
        failExit(std_ext::formatter() << "Cannot find tree "
                 << treeName << " in file " << cmd_input->getValue());
    tree.LinkBranches(tree.Tree);


    histOut = cmd_output->isSet();
    unique_ptr<WrapTFileOutput> outFile;
    if(histOut) {
        outFile = std_ext::make_unique<WrapTFileOutput>(
                      cmd_output->getValue(),
                      WrapTFileOutput::mode_t::recreate,
                      true);
    }


    HistogramFactory histfac("TriplePi0_plot");

    auto makeProbHist = [&histfac](const std::string& title)
    {
        BinSettings bs(200,0,1);
        return histfac.makeTH2D(title, "prob_{sig}","prob_{bkg}",bs,bs);
    };
    auto getAfterProbCut = [&histfac](TH2D* probs)
    {
        auto bsx = getBins(probs->GetXaxis());
        auto bsy = getBins(probs->GetYaxis());
        auto retH = histfac.makeTH2D(std_ext::formatter() << "sig and anit-bkg: " << probs->GetTitle(),
                                     "prob_{sig}","prob_{bkg}",
                                     bsx,bsy);
        auto maxx = probs->GetNbinsX();
        auto maxy = probs->GetNbinsY();

        for (int binx = 1 ; binx <= maxx ; ++binx)
        {
            for (int biny = 1 ; biny <= maxy; ++biny)
            {
                double acc = 0;
                for ( int ix = binx; ix <= maxx ; ++ix)
                    for ( int iy = 1 ; iy < biny ; ++iy)
                        acc += probs->GetBinContent(ix,iy);
                retH->SetBinContent(binx,biny,acc);
            }
        }
        return retH;
    };
    auto hProbs    = makeProbHist("All events");
    auto hProbsSig = makeProbHist("Signal");
    auto hProbsBkg = makeProbHist("Background");



    auto makeProtonHist = [&histfac](const string& title)
    {
        const string eKinProton("E^{kin}_{proton} [MeV]");
        const string thetaProton("#theta^{kin}_{proton} [#circ]");
        return histfac.makeTH2D(title,eKinProton,thetaProton,BinSettings(350,0,100),BinSettings(300,0,75));
    };
    auto hProtonSelection = makeProtonHist("kin-fitted Proton");

    for ( auto en=0u ; en < tree.Tree->GetEntries() ; ++en)
    {
        tree.Tree->GetEntry(en);


        hProbs->Fill(tree.SIG_prob,tree.BKG_prob);
        if (tree.MCTrue == 1)
            hProbsSig->Fill(tree.SIG_prob,tree.BKG_prob);
        if (tree.MCTrue == 2)
            hProbsBkg->Fill(tree.SIG_prob,tree.BKG_prob);

        hProtonSelection->Fill(tree.EMB_proton().E() - 938.3,std_ext::radian_to_degree(tree.EMB_proton().Theta()));
    }

    auto hAfterProbCuts    = getAfterProbCut(hProbs);
    auto hAfterProbCutsSig = getAfterProbCut(hProbsSig);
    auto hAfterProbCutsBkg = getAfterProbCut(hProbsBkg);



    // OUTPUT ==============================

    argc=1; // prevent TRint to parse any cmdline except prog name
    auto app = cmd_batchmode->isSet() || !std_ext::system::isInteractive()
               ? nullptr
               : std_ext::make_unique<TRint>("Ant-makeSigmas",&argc,argv,nullptr,0,true);
    if(app) {
        auto colz = drawoption("colz");

        canvas("probs and cuts")
                << colz
                << hProbs << hProbsSig << hProbsBkg
                << hAfterProbCuts << hAfterProbCutsSig << hAfterProbCutsBkg
                << endc;

        canvas("proton")
                << colz
                << hProtonSelection
                << endc;


        if(outFile)
            LOG(INFO) << "Stopped running, but close ROOT properly to write data to disk.";
        app->Run(kTRUE); // really important to return...
        if(outFile)
            LOG(INFO) << "Writing output file...";
        outFile = nullptr;   // and to destroy the master WrapTFile before TRint is destroyed
    }

    return EXIT_SUCCESS;
}



