#include <iostream>
#include <map>
#include <list>
#include <memory>

#include <base/std_ext.h>
#include "calibration/gui/FitCanvas.h"
#include "calibration/gui/FitFunction.h"
#include "base/Logger.h"
#include "base/CmdLine.h"

#include "calibration/gui/FitGausPol3.h"

#include "TRint.h"
#include "TH1D.h"
#include "TFile.h"

using namespace std;
using namespace ant;
using namespace ant::calibration::gui;


int main(int argc, char** argv) {
    SetupLogger();

    TCLAP::CmdLine cmd("canvas_test", ' ', "0.1");

    auto filename = cmd.add<TCLAP::ValueArg<string>>("f","file","ROOT File to open",true,"","inputfile");
    auto histname = cmd.add<TCLAP::ValueArg<string>>("H","hist","ROOT Histogram name in file",true,"","histogram");

    cmd.parse(argc, argv);

//    if(!filename->isSet() || !histname->isSet())
//        return 1;

    el::Loggers::setVerboseLevel(9);

    int i=0;
    TRint app("omega",&i,nullptr);

    CalCanvas* c = new CalCanvas("test");

    TFile* f = new TFile(filename->getValue().c_str(),"READ");
    TH1D* h = nullptr;
    f->GetObject(histname->getValue().c_str(), h);
    auto f1 = std::make_shared<FitGausPol3>();
    c->Show(h,f1.get());
    app.Run(kFALSE);

}
