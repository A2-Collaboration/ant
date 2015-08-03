#include <iostream>
#include <map>
#include <list>
#include <memory>

#include <base/std_ext.h>
#include "calibration/gui/FitCanvas.h"
#include "calibration/gui/FitFunction.h"
#include "base/Logger.h"

#include "calibration/gui/FitGausPol3.h"

#include "TRint.h"
#include "TH1D.h"

using namespace std;
using namespace ant::calibration::gui;


int main(int argc, char** argv) {
    SetupLogger();
    el::Loggers::setVerboseLevel(9);

    TRint app("omega",&argc,argv);

    CalCanvas* c = new CalCanvas("test");
    TH1D* h = new TH1D("h","h",100,-10,10);
    h->FillRandom("gaus",1000);

    auto f1 = std::make_shared<FitGausPol3>();
    c->Show(h,f1);
    app.Run(kTRUE);

}

