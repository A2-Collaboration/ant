#include "detail/TPCSim_tools.h"

#include "base/Logger.h"
#include "tclap/CmdLine.h"

#include "base/std_ext/math.h"

#include "TSystem.h"
#include "TRint.h"

#include "TGraphErrors.h"
#include "TGraph.h"
#include "TApplication.h"
#include "TAxis.h"
#include "TH2D.h"

#include <iostream>

using namespace std;
using namespace ant;
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

    cmd.parse(argc, argv);
    if(cmd_verbose->isSet()) {
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());
    }

    const TPCSim::resolution_t single_point_res = {0.075,0.053}; // (x,z)
    const TPCSim::tpcproperties tpc;

    const auto p = TPCSim::generatePoints( 0.0, std_ext::degree_to_radian(20.0), single_point_res, tpc);




    LOG(INFO) << p;

    argc=0; // prevent TRint to parse any cmdline
    TRint app(argv[0], &argc, argv, nullptr, 0, true);


    auto canvas = new TH2D("","", 1,-30, 30, 1, -30 ,30);
    canvas->SetStats(false);
    canvas->Draw();

    auto tpcarea = tpc.getOutline();
    tpcarea->Draw("L same");
    tpcarea->GetXaxis()->SetTitle("r [cm]");
    tpcarea->GetXaxis()->SetRangeUser(-3,15);
    tpcarea->GetYaxis()->SetTitle("z [cm]");

    auto g = TPCSim::makeGraph(p,single_point_res,tpc);
    g->Draw("P same");

    auto target = [] () {
        auto g = new TGraph(5);
        constexpr auto l=10.0;
        constexpr auto d=2.0;
        g->SetPoint(0,d/2,l/2);
        g->SetPoint(1,-d/2,l/2);
        g->SetPoint(2,-d/2,-l/2);
        g->SetPoint(3,d/2,-l/2);
        g->SetPoint(4,d/2,l/2);
        return g;
    }();
    target->Draw("L same");

    auto cb = [] () {
        constexpr auto np = 180;
        constexpr auto r = 24.5;
        auto g = new TGraph(np+1);
        for(int i=0; i<np; ++i) {
            const auto phi = std_ext::degree_to_radian(360.0/np*i);
            g->SetPoint(i,r*cos(phi),r*sin(phi));
        }
        g->SetPoint(np, r, 0.0);
        return g;
    }();
    cb->Draw("L same");

    ///@todo: draw fitted result as line

    app.Run(kTRUE); // really important to return...


    return EXIT_SUCCESS;
}
