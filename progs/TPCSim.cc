#include "detail/TPCSim_tools.h"

#include "base/Logger.h"
#include "tclap/CmdLine.h"

#include "base/std_ext/math.h"

#include "TSystem.h"
#include "TRint.h"

#include "TGraphErrors.h"
#include "TGraph.h"
#include "TF1.h"
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

    const auto points = TPCSim::generatePoints( 0.0, std_ext::degree_to_radian(45.0), single_point_res, tpc);

    auto makeFit = [](const vector<vec2>& points, const resolution_t& res, const tpcproperties& tpc)
    {
    vector<Value_t> rs;
    vector<Value_t> zs;

    for (const auto& point: points)
    {
        const auto uncert = getUncertainties(point, res, tpc);
        rs.push_back(Value_t{point.x,uncert.x});
        zs.push_back(Value_t{point.y,uncert.y});
    }

    return trackFitter_t(rs,zs);
    };

    auto fit = makeFit(points,single_point_res,tpc);

    LOG(INFO) << points;
    LOG(INFO) << "fitted rs: " << fit.Fitted_Rs;
    LOG(INFO) << "fitted zs: " << fit.Fitted_Zs;


    argc=0; // prevent TRint to parse any cmdline
    TRint app(argv[0], &argc, argv, nullptr, 0, true);


    TPCSim::draw::makeCanvas(tpc)->Draw();

    auto scene = TPCSim::draw::makeScene(tpc);
    for_each(scene.begin(),scene.end(),[](TGraph* g){g->Draw("L same");});

    auto g = TPCSim::draw::makeGraph(points,single_point_res,tpc);
    g->Draw("P same");

    auto fitfkt = TPCSim::draw::makeFitTF1(fit);
    fitfkt->Draw("same");


    app.Run(kTRUE); // really important to return...


    return EXIT_SUCCESS;
}
