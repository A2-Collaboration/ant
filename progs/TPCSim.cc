#include "detail/TPCSim_tools.h"

#include "base/Logger.h"
#include "tclap/CmdLine.h"

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
    auto histSigmaZ0    = histfac.makeTH2D("","#theta [#circ]","#Delta_{z(0)} [cm]",thetabins,deltsZbins,"sigmaZ0");
    auto histSigmaTheta = histfac.makeTH2D("","#theta [#circ]","#Delta_{#theta} [#circ]",thetabins,deltaThetabins,"sigmaTheta");
    trackFitter_t fitter;

    const unsigned nTracks  = 100000;
    const unsigned nHitCut  = 3;
    using trackHits_t = vector<ant::vec2>;
    vector<trackHits_t> simTracks;
    vector<trackFitter_t::result_t> fitResults;
    simTracks.reserve(nTracks);
    fitResults.reserve(nTracks);

    TRandom3 rnd;
    rnd.SetSeed(0);

    for (auto ntrack = 0u; ntrack < nTracks; ++ntrack)
    {
        const auto angle = rnd.Uniform(M_PI);
        const auto z0    = rnd.Uniform(10.) - 5.;
        const auto points = TPCSim::generatePoints( z0, angle, single_point_res, tpc);
        if(points.size() < nHitCut) {continue;}
        const auto result = fitter.DoFit(points,single_point_res,tpc);
        histSigmaZ0->Fill(radian_to_degree(angle),result.GetDZ0(z0));
        histSigmaTheta->Fill(radian_to_degree(angle),radian_to_degree(result.GetDTheta(angle)));

        if (show_debug)
        {
            simTracks.push_back(points);
            fitResults.push_back(result);
        }
    }

    argc=0; // prevent TRint to parse any cmdline
    TRint app(argv[0], &argc, argv, nullptr, 0, true);

    if (show_debug)
    {
        TPCSim::draw::makeCanvas(tpc)->Draw();

        auto scene = TPCSim::draw::makeScene(tpc);
        for_each(scene.begin(),scene.end(),[](TGraph* g){g->Draw("L same");});

        for_each(simTracks.begin(),simTracks.end(),[&single_point_res,&tpc](const trackHits_t& tr){draw::makeGraph(tr,single_point_res,tpc)->Draw("P same");});


        for_each(fitResults.begin(),fitResults.end(),[](const trackFitter_t::result_t& res){draw::makeFitTF1(res)->Draw("same");});
    }


    new TCanvas();
    histSigmaTheta->Draw("colz");

    new TCanvas();
    histSigmaZ0->Draw("colz");

    app.Run(kTRUE); // really important to return...


    return EXIT_SUCCESS;
}
