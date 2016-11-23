#include <string>
#include <map>

#include "analysis/physics/production/triplePi0.h"
#include "analysis/plot/root_draw.h"

#include "base/CmdLine.h"
#include "base/Logger.h"
#include "base/PlotExt.h"
#include "base/std_ext/math.h"
#include "base/std_ext/string.h"
#include "base/std_ext/system.h"
#include "base/WrapTFile.h"

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



static volatile bool interrupt = false;
static bool noStore = false;
static bool histOut = false;


auto failExit = [] (const string& message)
{
    LOG(ERROR) << message;
    exit(EXIT_FAILURE);
};


struct sigOverBkg_t
{
    shared_ptr<HistogramFactory> HistFac;
    BinSettings Bs;
    TH2D* probs;
    TH2D* SigC;
    TH2D* Ratio() const
    {

        auto retH = HistFac->makeTH2D("sig/bkg","prob_{sig}","prob_{bkg}",Bs,Bs);
        return retH;
    }
    sigOverBkg_t(const shared_ptr<HistogramFactory>& histFac, const BinSettings& bs):
        HistFac(histFac),
        Bs(bs)
    {
        SigC = HistFac->makeTH2D("signal","prob_{sig}","prob_{bkg}",Bs,Bs);
        probs = HistFac->makeTH2D("background","prob_{sig}","prob_{bkg}",Bs,Bs);
    }

    void fillProbs(double pSig, double pBkg)
    {
        probs->Fill(pSig,pBkg);
    }
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
    auto cmd_output     = cmd.add<TCLAP::ValueArg<string>>("",  "save-hist", "Save results to a histogram", false, "","filename");

    //switches
    auto cmd_batchmode  = cmd.add<TCLAP::SwitchArg>("b", "batch",     "Run in batch mode (no ROOT shell afterwards)");

    cmd.parse(argc, argv);


    WrapTFileInput inFile(cmd_input->getValue());
    auto treeName = triplePi0::treeAccessName();
    triplePi0::PionProdTree tree;
    if (!inFile.GetObject(treeName,tree.Tree))
        failExit(std_ext::formatter() << "Cannot find tree " << treeName << " in file " << cmd_input->getValue());
    tree.LinkBranches(tree.Tree);


    histOut = cmd_output->isSet();
    unique_ptr<WrapTFileOutput> outFile;
    if(histOut) {
        outFile = std_ext::make_unique<WrapTFileOutput>(cmd_output->getValue(),
                                                           WrapTFileOutput::mode_t::recreate,
                                                           true);
    }


    auto histfac = make_shared<HistogramFactory>("TriplePi0_plot");

    BinSettings bsProb(100,0,1);

    sigOverBkg_t sigBkg(histfac,bsProb);


    auto makeProtonHist = [&histfac](const string& title)
    {
        const string eKinProton("E^{kin}_{proton} [MeV]");
        const string thetaProton("#theta^{kin}_{proton} [#circ]");
        return histfac->makeTH2D(title,eKinProton,thetaProton,BinSettings(350,0,100),BinSettings(300,0,75));
    };
    auto hProtonSelection = makeProtonHist("kin-fitted Proton");

    for ( auto en=0u ; en < tree.Tree->GetEntries() ; ++en)
    {
        tree.Tree->GetEntry(en);

        sigBkg.fillProbs(tree.SIG_prob,tree.BKG_prob);

        hProtonSelection->Fill(tree.EMB_proton().E() - 938.3,std_ext::radian_to_degree(tree.EMB_proton().Theta()));
    }



    // OUTPUT ==============================

    argc=1; // prevent TRint to parse any cmdline except prog name
    auto app = cmd_batchmode->isSet() || !std_ext::system::isInteractive()
               ? nullptr
               : std_ext::make_unique<TRint>("Ant-makeSigmas",&argc,argv,nullptr,0,true);
    if(app) {



        if(outFile)
            LOG(INFO) << "Stopped running, but close ROOT properly to write data to disk.";
        app->Run(kTRUE); // really important to return...
        if(outFile)
            LOG(INFO) << "Writing output file...";
        outFile = nullptr;   // and to destroy the master WrapTFile before TRint is destroyed
    }

    return EXIT_SUCCESS;
}



