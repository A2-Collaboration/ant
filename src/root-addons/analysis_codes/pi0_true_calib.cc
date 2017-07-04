#include "pi0_true_calib.h"
#include "base/Logger.h"
#include "analysis/plot/HistogramFactory.h"

#include "TH1D.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;



void pi0_true_calib::Do()
{
    HistogramFactory HistFac("h");

    LOG(INFO) << "Hej!";
    auto h = HistFac.makeTH1D("bla",{"x",{100,-5,5}});
    h->Draw();

}
