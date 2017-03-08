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

const IntervalD protonIMcut(ParticleTypeDatabase::Proton.Mass() - 180,
                            ParticleTypeDatabase::Proton.Mass() + 180);
const double    ProbCut(0.1);




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
    auto nBins  = 100u;
    auto m3pi0  = histfac.makeTH1D("m(3#pi^0)","m(3#pi^0) [MeV]","#",      BinSettings(nBins,0,1000));
    auto mK0S   = histfac.makeTH1D("m(K^{0}_{S})","m(K^{0}_{S}) [MeV]","#",BinSettings(nBins,300,700));
    auto mpPi0  = histfac.makeTH1D("m(p #pi^{0})","m(p #pi^{0}) [MeV]","#",BinSettings(nBins,1100,1600));

    auto mK0SF  = histfac.makeTH1D("m(K^{0}_{S}) Fit","m(K^{0}_{S}) [MeV]","#",BinSettings(nBins,300,700));
    auto mpPi0F = histfac.makeTH1D("m(p #pi^{0}) Fit","m(p #pi^{0}) [MeV]","#",BinSettings(nBins,1100,1600));

    auto testCuts = [] (const triplePi0::PionProdTree& tree)
    {
        return (!protonIMcut.Contains(tree.proton_MM().M()) &&
                tree.EMB_prob  < 0.1 &&
                tree.SIG_prob  < 0.1 &&
                tree.IM6g      < 640.0  );
    };


    for ( auto en=0u ; en < tree.Tree->GetEntries() ; ++en)
    {
        tree.Tree->GetEntry(en);

//        if (testCuts(tree))
//                continue;

        if (tree.SIGMA_combination().size() == 6)
        {
            TLorentzVector vK0S(0,0,0,0);
            TLorentzVector vpPi0(tree.proton());
            for (auto i = 0u; i < 4; ++i)
                vK0S += tree.photons().at(tree.SIGMA_combination().at(i));

            for (auto i = 4u; i < 6; ++i)
                vpPi0 += tree.photons().at(tree.SIGMA_combination().at(i));
            mK0S->Fill(vK0S.M());
            mpPi0->Fill(vpPi0.M());
        }

        mK0SF->Fill(tree.SIGMA_k0s().M());
        mpPi0F->Fill(tree.SIGMA_SigmaPlus().M());

        m3pi0->Fill(tree.IM6g);




    }


    // OUTPUT ==============================

    argc=1; // prevent TRint to parse any cmdline except prog name
    auto app = cmd_batchmode->isSet() || !std_ext::system::isInteractive()
               ? nullptr
               : std_ext::make_unique<TRint>("Ant-makeSigmas",&argc,argv,nullptr,0,true);
    if(app) {
        auto colz = drawoption("colz");

        canvas("hists")
                << m3pi0
                << mK0S
                << mpPi0
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



