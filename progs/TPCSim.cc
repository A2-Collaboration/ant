#include "detail/TPCSim_tools.h"

#include "base/Logger.h"
#include "tclap/CmdLine.h"

#include "base/std_ext/math.h"

#include "TSystem.h"
#include "TRint.h"

#include <iostream>

using namespace std;
using namespace ant;
using namespace ant::TPCSim;

static volatile bool interrupt = false;

auto residuals = [] (const Value_t& a, const Value_t& b, const vector<Value_t>& z, const vector<Value_t>& r) {
    vector<double> residuals(r.size());
    transform(z.begin(), z.end(), r.begin(), residuals.begin(),
              [&a, &b] (const double& x_i, const double& y_i) {
        return a + b*x_i - y_i;
    });
    return residuals;
};

int main(const int argc, const char** argv) {
    SetupLogger();

    signal(SIGINT, [] (int) { interrupt = true;} );

    TCLAP::CmdLine cmd("TPCSim", ' ', "1");
    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"int");

    cmd.parse(argc, argv);
    if(cmd_verbose->isSet()) {
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());
    }

    const TPCSim::driftproperties d = {400,400};
    const TPCSim::tpcproperties tpc;

    TPCSim::generatePoints(0,std_ext::degree_to_radian(45.0), d, tpc);


    return EXIT_SUCCESS;
}
