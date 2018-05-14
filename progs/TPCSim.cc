#include "detail/TPCSim_tools.h"

#include "base/Logger.h"
#include "tclap/CmdLine.h"
#include "base/ForLoopCounter.h"

#include "base/std_ext/math.h"
#include "analysis/plot/HistogramFactory.h"


#include "TSystem.h"
#include "TRint.h"

#include "TGraphErrors.h"
#include "TGraph.h"
#include "TF1.h"
#include "TApplication.h"
#include "TAxis.h"
#include "TH2D.h"
#include "TRandom3.h"
#include "TCanvas.h"

#include <iostream>

using namespace std;
using namespace ant;
using namespace ant::std_ext;
using namespace TPCSim;

static volatile bool interrupt = false;



ostream& operator<<(ostream& o, const vec2& v) {
  o << "(" << v.x <<"," << v.y << ")";
  return o;
}


int main(int argc, char** argv) {
    SetupLogger();

    signal(SIGINT, [] (int) { interrupt = true;} );

    TCLAP::CmdLine cmd("TPCSim", ' ', "1");
    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"int");
    auto cmd_debug_plots = cmd.add<TCLAP::MultiSwitchArg>("d","debug-plot","show scatch of all generated hits and fitfunctions",false);


    cmd.parse(argc, argv);
    if(cmd_verbose->isSet()) {
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());
    }

    const auto show_debug = cmd_debug_plots->getValue();

    const TPCSim::resolution_t single_point_res = {0.075,0.053}; // (x,z)
    const TPCSim::tpcproperties tpc;

    ant::analysis::HistogramFactory histfac("tpcsim");
    const ant::BinSettings thetabins(100,0,180);
    const ant::BinSettings deltsZbins(100,-1,1);
    const ant::BinSettings deltaThetabins(100,-5,5);
    // quick check:
    auto histSigmaZ0    = histfac.makeTH2D("","#theta [#circ]","#Delta_{z(0)} [cm]",thetabins,deltsZbins,"sigmaZ0");
    auto histSigmaTheta = histfac.makeTH2D("","#theta [#circ]","#Delta_{#theta} [#circ]",thetabins,deltaThetabins,"sigmaTheta");

    tree_t tree;
    tree.CreateBranches(histfac.makeTTree(tree_t::treeName()));


    trackFitter_t fitter;

    const unsigned nTracks  = 100000;
//    const unsigned nHitCut  = 3;
    using trackHits_t = vector<ant::vec2>;
    vector<trackHits_t> simTracks;
    vector<trackFitter_t::result_t> fitResults;
    simTracks.reserve(nTracks);
    fitResults.reserve(nTracks);

    TRandom3 rnd;
    rnd.SetSeed(0);

    for (auto nHitCut = 2; nHitCut <= tpc.nRings ; ++nHitCut)
    {
        LOG(INFO) << "Hits per Track >= " << nHitCut;
        for (const auto track: ForLoopCounter<size_t>(nTracks))
        {
            const auto angle = rnd.Uniform(M_PI);
            const auto z0    = rnd.Uniform(10.) - 5.;
            const auto points = TPCSim::generatePoints( z0, angle, single_point_res, tpc);

            if((int)points.size() < nHitCut) {continue;}
            const auto result = fitter.DoFit(points,single_point_res,tpc);

            const auto dz0 = result.GetDZ0(z0);
            const auto dtheta = radian_to_degree(result.GetDTheta(angle));
            const auto angleDeg = radian_to_degree(angle);

            histSigmaZ0->Fill(angleDeg,dz0);
            histSigmaTheta->Fill(angleDeg,dtheta);

            tree.dTheta = dtheta;
            tree.dZ0 = dz0;
            tree.genTheta = angleDeg;
            tree.genZ0 = z0;
            tree.nHits = points.size();
            tree.nHitsCut = nHitCut;

            tree.Tree->Fill();

            if (show_debug)
            {
                simTracks.push_back(points);
                fitResults.push_back(result);
            }
        }
    }

    argc=0; // prevent TRint to parse any cmdline
    TRint app(argv[0], &argc, argv, nullptr, 0, true);

    if (show_debug)
    {
        TPCSim::draw::makeCanvas(tpc)->Draw();

        auto scene = TPCSim::draw::makeScene(tpc);
        for_each(scene.begin(),scene.end(),[](TGraph* g){g->Draw("L same");});

        for_each(simTracks.begin(),simTracks.end(),
                 [&single_point_res,&tpc](const trackHits_t& tr){ draw::makeGraph(tr,single_point_res,tpc)->Draw("P same"); });

        for_each(fitResults.begin(),fitResults.end(),
                 [](const trackFitter_t::result_t& res){ draw::makeFitTF1(res)->Draw("same"); });
    }


    auto c_theta = new TCanvas();
    c_theta->cd();
    histSigmaTheta->Draw("colz");

    auto c_sigma = new TCanvas();
    c_sigma->cd();
    histSigmaZ0->Draw("colz");

    app.Run(kTRUE); // really important to return...


    return EXIT_SUCCESS;
}
